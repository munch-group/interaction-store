"""
indra_cache.py
--------------
Persistent cache for INDRA DB REST API queries.

Stores results per-gene in indra_db_cache.json so repeated queries
are instant. Both the notebook and MCP server import from here.

Cache structure:
{
  "GENE_NAME": {
    "fetched": "2026-03-12",
    "ev_limit": 5,
    "statements": [ ... INDRA JSON ... ]
  },
  ...
}

Usage:
    from indra_cache import cached_get_statements, load_cache

    # Single gene (checks cache first, queries API on miss)
    stmts = cached_get_statements('ADRA2C', ev_limit=5)

    # Batch
    all_stmts = cached_get_statements_batch(['ADRA2C', 'IRS2', 'MAPT'])
"""

import json
import pathlib
import time
from datetime import date
from typing import Optional

CACHE_PATH = pathlib.Path(__file__).parent / 'indra_db_cache.json'


def load_cache() -> dict:
    if not CACHE_PATH.exists():
        return {}
    with open(CACHE_PATH) as f:
        data = json.load(f)
    if not isinstance(data, dict):
        # Cache file was corrupted (e.g. overwritten with a flat list).
        # Move it aside and start fresh so we don't block the user.
        backup = CACHE_PATH.with_suffix('.json.bak')
        CACHE_PATH.rename(backup)
        print(f'WARNING: indra_db_cache.json was not a dict — '
              f'moved to {backup.name} and starting fresh cache.')
        return {}
    return data


def save_cache(cache: dict):
    with open(CACHE_PATH, 'w') as f:
        json.dump(cache, f, indent=2)


def cache_age_days(gene: str, cache: Optional[dict] = None) -> Optional[int]:
    """Return age of cached entry in days, or None if not cached."""
    cache = cache or load_cache()
    entry = cache.get(gene)
    if not entry:
        return None
    from datetime import datetime
    fetched = datetime.strptime(entry['fetched'], '%Y-%m-%d').date()
    return (date.today() - fetched).days


def cached_get_statements(
    gene: str,
    ev_limit: int = 5,
    max_age_days: Optional[int] = None,
) -> list:
    """
    Get INDRA statements for a gene, using cache when available.

    Parameters
    ----------
    gene : str
        Gene name to query.
    ev_limit : int
        Max evidence objects per statement (passed to INDRA API).
    max_age_days : int, optional
        If set, re-query if cached result is older than this many days.
        If None, cached results never expire (delete manually to refresh).

    Returns
    -------
    list of INDRA Statement objects
    """
    from indra.statements import stmts_from_json, stmts_to_json

    cache = load_cache()
    entry = cache.get(gene)

    # Check cache hit
    if entry is not None:
        expired = False
        if max_age_days is not None:
            age = cache_age_days(gene, cache)
            expired = age is not None and age > max_age_days
        if not expired:
            return stmts_from_json(entry['statements'])

    # Cache miss — query API
    from indra.sources.indra_db_rest import get_statements
    proc = get_statements(agents=[gene], ev_limit=ev_limit)
    stmts = proc.statements

    # Store in cache
    cache[gene] = {
        'fetched': str(date.today()),
        'ev_limit': ev_limit,
        'statements': stmts_to_json(stmts),
    }
    save_cache(cache)

    return stmts


def cached_get_statements_batch(
    genes: list[str],
    ev_limit: int = 5,
    max_age_days: Optional[int] = None,
    sleep_between: float = 1.0,
    verbose: bool = True,
) -> list:
    """
    Query multiple genes, using cache where possible.

    Returns combined list of all statements. Prints progress if verbose.
    """
    all_stmts = []
    cache = load_cache()

    for gene in genes:
        entry = cache.get(gene)
        hit = False
        if entry is not None:
            expired = False
            if max_age_days is not None:
                age = cache_age_days(gene, cache)
                expired = age is not None and age > max_age_days
            if not expired:
                hit = True

        if hit:
            from indra.statements import stmts_from_json
            stmts = stmts_from_json(entry['statements'])
            if verbose:
                print(f'  {gene}: {len(stmts)} statements (cached)')
        else:
            try:
                stmts = cached_get_statements(gene, ev_limit=ev_limit,
                                              max_age_days=max_age_days)
                if verbose:
                    print(f'  {gene}: {len(stmts)} statements (fetched)')
                time.sleep(sleep_between)
            except Exception as e:
                if verbose:
                    print(f'  {gene}: FAILED ({e})')
                stmts = []

        all_stmts.extend(stmts)

    if verbose:
        print(f'\nTotal: {len(all_stmts)} statements '
              f'({sum(1 for g in genes if g in load_cache())} cached)')

    return all_stmts


def cache_summary() -> str:
    """Return a summary of what's in the cache."""
    cache = load_cache()
    if not cache:
        return 'Cache is empty.'
    lines = [f'INDRA DB cache: {len(cache)} genes\n']
    total = 0
    for gene in sorted(cache):
        entry = cache[gene]
        n = len(entry.get('statements', []))
        total += n
        age = cache_age_days(gene, cache)
        lines.append(f'  {gene:12s}  {n:5d} stmts  fetched {entry["fetched"]}'
                      f'  ({age}d ago)')
    lines.append(f'\nTotal: {total} cached statements')
    return '\n'.join(lines)


def invalidate_cache(genes: Optional[list[str]] = None):
    """Remove genes from cache. If genes is None, clear everything."""
    if genes is None:
        if CACHE_PATH.exists():
            CACHE_PATH.unlink()
        return
    cache = load_cache()
    for g in genes:
        cache.pop(g, None)
    save_cache(cache)
