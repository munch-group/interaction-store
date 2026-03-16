"""
indra_cache.py
--------------
Persistent cache for INDRA DB REST API queries.

Stores results per-gene in a SQLite database (indra_db_cache.db) so
repeated queries are instant. Both the notebook and MCP server import
from here.

Schema:
    CREATE TABLE cache (
        gene       TEXT PRIMARY KEY,
        fetched    TEXT NOT NULL,        -- ISO date YYYY-MM-DD
        ev_limit   INTEGER NOT NULL,
        statements TEXT NOT NULL         -- JSON array (INDRA statement list)
    );

Usage:
    from indra_cache import cached_get_statements

    # Single gene (checks cache first, queries API on miss)
    stmts = cached_get_statements('ADRA2C', ev_limit=5)

    # Multiple genes
    all_stmts = cached_get_statements(['ADRA2C', 'IRS2', 'MAPT'], verbose=True)
"""

import json
import pathlib
import sqlite3
import time
from datetime import date, datetime
from typing import Optional
from tqdm.auto import tqdm

CACHE_JSON_PATH = pathlib.Path(__file__).parent / 'indra_db_cache.json'
CACHE_DB = pathlib.Path(__file__).parent / 'indra_db_cache.db'

# Keep old path reference for backward compat
CACHE_PATH = CACHE_JSON_PATH


# ---------------------------------------------------------------------------
# SQLite helpers
# ---------------------------------------------------------------------------

def _get_db() -> sqlite3.Connection:
    """Return a sqlite3 connection, creating table if needed."""
    conn = sqlite3.connect(CACHE_DB)
    conn.execute('''CREATE TABLE IF NOT EXISTS cache (
        gene       TEXT PRIMARY KEY,
        fetched    TEXT NOT NULL,
        ev_limit   INTEGER NOT NULL,
        statements TEXT NOT NULL
    )''')
    return conn


def _get_entry(gene: str) -> dict | None:
    """Fetch a single gene's cache entry."""
    conn = _get_db()
    row = conn.execute(
        'SELECT fetched, ev_limit, statements FROM cache WHERE gene = ?',
        (gene,)
    ).fetchone()
    conn.close()
    if row is None:
        return None
    return {'fetched': row[0], 'ev_limit': row[1], 'statements': json.loads(row[2])}


def _put_entry(gene: str, fetched: str, ev_limit: int, statements_json: list):
    """Upsert a single gene's cache entry."""
    conn = _get_db()
    conn.execute(
        'INSERT OR REPLACE INTO cache (gene, fetched, ev_limit, statements) VALUES (?, ?, ?, ?)',
        (gene, fetched, ev_limit, json.dumps(statements_json))
    )
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------------
# Migration: JSON → SQLite (one-time, runs at import)
# ---------------------------------------------------------------------------

def _migrate_json_to_db():
    """Convert indra_db_cache.json to indra_db_cache.db (one-time)."""
    if not CACHE_JSON_PATH.exists() or CACHE_DB.exists():
        return
    print(f'Migrating {CACHE_JSON_PATH.name} → {CACHE_DB.name} ...')
    try:
        with open(CACHE_JSON_PATH) as f:
            data = json.load(f)
    except (json.JSONDecodeError, ValueError) as exc:
        print(f'WARNING: {CACHE_JSON_PATH.name} is corrupt ({exc.__class__.__name__}) '
              f'— skipping migration.')
        return
    if not isinstance(data, dict):
        print(f'WARNING: {CACHE_JSON_PATH.name} was not a dict — skipping migration.')
        return
    conn = _get_db()
    for gene, entry in data.items():
        conn.execute(
            'INSERT OR REPLACE INTO cache (gene, fetched, ev_limit, statements) VALUES (?, ?, ?, ?)',
            (gene, entry['fetched'], entry['ev_limit'], json.dumps(entry['statements']))
        )
    conn.commit()
    conn.close()
    print(f'Migrated {len(data)} genes. You can now delete {CACHE_JSON_PATH.name}.')


_migrate_json_to_db()


# ---------------------------------------------------------------------------
# Public API (signatures unchanged)
# ---------------------------------------------------------------------------

def load_cache() -> dict:
    """Load entire cache as dict. Prefer _get_entry() for single-gene access."""
    conn = _get_db()
    rows = conn.execute('SELECT gene, fetched, ev_limit, statements FROM cache').fetchall()
    conn.close()
    return {
        row[0]: {'fetched': row[1], 'ev_limit': row[2], 'statements': json.loads(row[3])}
        for row in rows
    }


def cache_age_days(gene: str) -> Optional[int]:
    """Return age of cached entry in days, or None if not cached."""
    entry = _get_entry(gene)
    if not entry:
        return None
    fetched = datetime.strptime(entry['fetched'], '%Y-%m-%d').date()
    return (date.today() - fetched).days


_CONFIDENCE_LEVELS = {
    'curated':  'Database-curated interactions only (Signor, BioGRID, Reactome, etc.)',
    'balanced': 'Database-curated OR NLP-extracted with >= 5 independent evidence texts',
    'supported': 'NLP-extracted with >= 5 independent evidence texts',
}


def _build_query(gene: str, confidence: Optional[str],
                 namespace: Optional[str] = None):
    """Return an INDRA DB REST query object for the given gene and confidence.

    Returns None for the default (no confidence filter and no namespace) —
    caller should use the simple ``get_statements(agents=[g])`` path instead.
    """
    from indra.sources.indra_db_rest.query import (
        HasAgent, HasDatabases, HasEvidenceBound, Union,
    )

    has_filter = confidence is not None or namespace is not None
    if not has_filter:
        return None

    agent = HasAgent(gene, namespace=namespace) if namespace else HasAgent(gene)

    if confidence is None:
        return agent
    elif confidence == 'curated':
        return agent & HasDatabases()
    elif confidence == 'balanced':
        return agent & Union([HasDatabases(), HasEvidenceBound(['>= 5'])])
    elif confidence == 'supported':
        return agent & HasEvidenceBound(['>= 5'])
    else:
        raise ValueError(
            f"confidence must be one of {list(_CONFIDENCE_LEVELS)} or None, "
            f"got {confidence!r}"
        )


def _cache_key(gene: str, confidence: Optional[str],
               namespace: Optional[str] = None) -> str:
    """Return the cache key for a gene + confidence combination."""
    parts = [gene]
    if confidence:
        parts.append(confidence)
    if namespace:
        parts.append(f'ns={namespace}')
    return '::'.join(parts)


def _fetch_one(gene: str, confidence, ev_limit: int,
               namespace: Optional[str] = None) -> list:
    """Fetch statements for a single gene from the INDRA REST API.

    Returns a list of INDRA Statement objects (may be empty on failure).
    Raises on network/API errors so the caller can handle retries.

    Stderr is suppressed during the call because INDRA's internal
    pagination thread prints raw tracebacks on connection resets that
    cannot be caught from Python.
    """
    import sys, io
    query = _build_query(gene, confidence, namespace=namespace)
    # Suppress INDRA's internal thread tracebacks on connection resets.
    _real_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        if query is None:
            from indra.sources.indra_db_rest import get_statements
            proc = get_statements(agents=[gene], ev_limit=ev_limit)
            return proc.statements
        else:
            from indra.sources.indra_db_rest.api import get_statements_from_query
            proc = get_statements_from_query(query, ev_limit=ev_limit)
            return proc.statements if hasattr(proc, 'statements') else proc
    finally:
        sys.stderr = _real_stderr


def _fetch_with_retry(gene: str, confidence, ev_limit: int,
                      namespace: Optional[str] = None,
                      max_retries: int = 3, base_delay: float = 5.0,
                      verbose: bool = False) -> 'tuple[str, list, str|None]':
    """Fetch with exponential backoff on transient failures.

    Returns (gene, statements, error_msg_or_None).
    """
    import requests

    for attempt in range(max_retries):
        try:
            stmts = _fetch_one(gene, confidence, ev_limit, namespace=namespace)
            return (gene, stmts, None)
        except (requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
                requests.exceptions.HTTPError,
                ConnectionResetError,
                OSError) as e:
            # Transient network / rate-limit error — back off and retry
            delay = base_delay * (2 ** attempt)
            if verbose:
                tqdm.write(f'  {gene}: attempt {attempt+1}/{max_retries} failed '
                           f'({type(e).__name__}), retrying in {delay:.0f}s...')
            time.sleep(delay)
        except Exception as e:
            # Non-transient error — don't retry
            return (gene, [], str(e))

    return (gene, [], f'Failed after {max_retries} retries')


def cached_get_statements(
    gene: 'str | list[str]',
    ev_limit: int = 5,
    max_age_days: Optional[int] = None,
    sleep_between: float = 1.0,
    verbose: bool = False,
    confidence: Optional[str] = None,
    parallel: Optional[int] = None,
    namespace: Optional[str] = None,
) -> list:
    """
    Get INDRA statements for one or more genes, using cache when available.

    Parameters
    ----------
    gene : str or list[str]
        Gene name(s) to query.
    ev_limit : int
        Max evidence objects per statement (passed to INDRA API).
    max_age_days : int, optional
        If set, re-query if cached result is older than this many days.
        If None, cached results never expire (delete manually to refresh).
    sleep_between : float
        Seconds to sleep between API calls (only for sequential queries).
    verbose : bool
        Print progress per gene.
    confidence : str, optional
        Filter by confidence level. Options:

        - ``None`` (default): no filter, return all statements.
        - ``'curated'``: database-curated interactions only (Signor,
          BioGRID, Reactome, PhosphoSitePlus, etc.).
        - ``'balanced'``: database-curated OR NLP-extracted with >= 5
          independent evidence texts.
        - ``'supported'``: NLP-extracted with >= 5 independent evidence
          texts (no database requirement).

        Results are cached separately per confidence level.
    namespace : str, optional
        Restrict agent grounding to a specific database namespace
        (e.g. ``'HGNC'`` for genes only, ``'CHEBI'`` for chemicals).
        Cached separately per namespace.
    parallel : int, optional
        Number of concurrent API requests for cache misses (default None
        = sequential). Uses threads with exponential backoff on transient
        failures (connection resets, timeouts, HTTP 429/503). Recommended
        value: 3-4. If the API throttles connections, concurrency is
        automatically reduced to 1 for remaining genes.

    Returns
    -------
    list of INDRA Statement objects
    """
    from indra.statements import stmts_from_json, stmts_to_json

    genes = [gene] if isinstance(gene, str) else list(gene)

    # Split into cache hits and misses
    cached_results = {}   # gene -> stmts
    needs_fetch = []      # genes that need API calls

    assert confidence in (None, 'curated', 'balanced', 'supported'), \
        f"confidence must be one of {list(_CONFIDENCE_LEVELS)} or None, got {confidence!r}"

    for g in genes:
        key = _cache_key(g, confidence, namespace)
        entry = _get_entry(key)
        hit = False
        if entry is not None:
            expired = False
            if max_age_days is not None:
                age = cache_age_days(key)
                expired = age is not None and age > max_age_days
            if not expired:
                hit = True
        if hit:
            cached_results[g] = stmts_from_json(entry['statements'])
            if verbose:
                print(f'  {g}: {len(cached_results[g])} statements (cached)')
        else:
            needs_fetch.append(g)

    # Fetch cache misses
    fetched_results = {}  # gene -> stmts

    if needs_fetch:
        def _do_fetch(g):
            """Fetch, cache, and return (gene, stmts, error)."""
            g, stmts, err = _fetch_with_retry(
                g, confidence, ev_limit, namespace=namespace,
                max_retries=3, verbose=verbose,
            )
            if not err:
                key = _cache_key(g, confidence, namespace)
                _put_entry(key, str(date.today()), ev_limit,
                           stmts_to_json(stmts))
            return (g, stmts, err)

        workers = min(parallel, len(needs_fetch)) if parallel and parallel > 1 else 1

        from tqdm.contrib.concurrent import thread_map
        results = thread_map(
            _do_fetch, needs_fetch,
            max_workers=workers,
            desc='Fetching from INDRA DB',
            # disable=not verbose,
            tqdm_class=tqdm,
        )

        for g, stmts, err in results:
            fetched_results[g] = [] if err else stmts

    # Reassemble in original gene order
    all_stmts = []
    for g in genes:
        all_stmts.extend(cached_results.get(g, fetched_results.get(g, [])))

    if verbose and len(genes) > 1:
        n_cached = len(cached_results)
        n_fetched = sum(1 for g in fetched_results if fetched_results[g])
        print(f'\nTotal: {len(all_stmts)} statements '
              f'({n_cached} cached, {n_fetched} fetched)')

    return all_stmts


# Backward compat alias
get_statements_batch = cached_get_statements


def cache_summary() -> str:
    """Return a summary of what's in the cache."""
    conn = _get_db()
    rows = conn.execute(
        'SELECT gene, fetched, ev_limit, json_array_length(statements) FROM cache ORDER BY gene'
    ).fetchall()
    conn.close()
    if not rows:
        return 'Cache is empty.'
    lines = [f'INDRA DB cache: {len(rows)} genes\n']
    total = 0
    for gene, fetched, ev_limit, n_stmts in rows:
        total += n_stmts
        fetched_date = datetime.strptime(fetched, '%Y-%m-%d').date()
        age = (date.today() - fetched_date).days
        lines.append(f'  {gene:12s}  {n_stmts:5d} stmts  fetched {fetched}'
                      f'  ({age}d ago)')
    lines.append(f'\nTotal: {total} cached statements')
    return '\n'.join(lines)


def invalidate_cache(genes: Optional[list[str]] = None):
    """Remove genes from cache. If genes is None, clear everything."""
    if genes is None:
        if CACHE_DB.exists():
            CACHE_DB.unlink()
        return
    conn = _get_db()
    conn.executemany('DELETE FROM cache WHERE gene = ?', [(g,) for g in genes])
    conn.commit()
    conn.close()
