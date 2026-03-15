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


def cached_get_statements(
    gene: 'str | list[str]',
    ev_limit: int = 5,
    max_age_days: Optional[int] = None,
    sleep_between: float = 1.0,
    verbose: bool = False,
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
        Seconds to sleep between API calls (only for multi-gene queries).
    verbose : bool
        Print progress per gene.

    Returns
    -------
    list of INDRA Statement objects
    """
    from indra.statements import stmts_from_json, stmts_to_json

    genes = [gene] if isinstance(gene, str) else list(gene)
    all_stmts = []
    n_cached = 0

    for g in tqdm(genes):
        entry = _get_entry(g)

        # Check cache hit
        hit = False
        if entry is not None:
            expired = False
            if max_age_days is not None:
                age = cache_age_days(g)
                expired = age is not None and age > max_age_days
            if not expired:
                hit = True

        if hit:
            stmts = stmts_from_json(entry['statements'])
            n_cached += 1
            if verbose:
                print(f'  {g}: {len(stmts)} statements (cached)')
        else:
            try:
                from indra.sources.indra_db_rest import get_statements
                proc = get_statements(agents=[g], ev_limit=ev_limit)
                stmts = proc.statements
                _put_entry(g, str(date.today()), ev_limit, stmts_to_json(stmts))
                n_cached += 1
                if verbose:
                    print(f'  {g}: {len(stmts)} statements (fetched)')
                if len(genes) > 1:
                    time.sleep(sleep_between)
            except Exception as e:
                if verbose:
                    print(f'  {g}: FAILED ({e})')
                stmts = []

        all_stmts.extend(stmts)

    if verbose and len(genes) > 1:
        print(f'\nTotal: {len(all_stmts)} statements ({n_cached} cached)')

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
