"""
statement_store.py
------------------
I/O and query functions for the INDRA statement store (statements.json).

Public API
----------
load_store(path=None)          → list[Statement]
save_store(stmts, path=None)   → None
add_statements(new, path=None) → list[Statement]
ev(text, ...)                  → Evidence
ag(name, ...)                  → Agent
summarize(path=None)           → None  (prints summary)
stmts_for(gene, path=None)     → list[Statement]
all_contexts(store_path=None)  → dict[str, int]
genes_by_context(context, ...) → GeneList
interactors(gene, ...)         → list[dict]
query_statements(**kwargs)     → list[dict]
"""

import json
import pathlib
from datetime import date
from typing import Any

from indra.statements import (
    Agent, Evidence,
    stmts_to_json, stmts_from_json,
)

STORE_PATH = pathlib.Path(__file__).parent / 'statements.json'


def _atomic_json_write(path: pathlib.Path, data, **kwargs):
    """Write JSON atomically: temp file + rename to avoid corruption."""
    import tempfile
    kwargs.setdefault('indent', 2)
    tmp_fd, tmp_path = tempfile.mkstemp(
        dir=path.parent, suffix='.tmp', prefix=f'.{path.stem}_')
    try:
        with open(tmp_fd, 'w') as f:
            json.dump(data, f, **kwargs)
        pathlib.Path(tmp_path).replace(path)
    except BaseException:
        pathlib.Path(tmp_path).unlink(missing_ok=True)
        raise


# ── I/O ───────────────────────────────────────────────────────────────────────

def load_store(path=None):
    """Load statements from JSON file, return list of INDRA Statement objects."""
    path = path or STORE_PATH
    if not path.exists():
        return []
    with open(path) as f:
        return stmts_from_json(json.load(f))


def save_store(stmts, path=None):
    """Serialise statement list to JSON file (atomic write)."""
    path = pathlib.Path(path or STORE_PATH)
    _atomic_json_write(path, stmts_to_json(stmts))
    print(f'Saved {len(stmts)} statements → {path}')


def add_statements(new_stmts, path=None):
    """Append new statements to the persistent store."""
    existing = load_store(path)
    combined = existing + new_stmts
    save_store(combined, path)
    return combined


# ── Convenience constructors ─────────────────────────────────────────────────

def ev(text, pmid=None, doi=None, context=None, hypothesis=False, direct=True):
    """Convenience: create a manual Evidence object.

    Args:
        text:       Mechanistic reasoning or supporting sentence.
        pmid:       PubMed ID string, e.g. '23640493'.
        doi:        DOI string, e.g. '10.1161/ATVBAHA.113.301522'.
        context:    Module name or thematic tag.
        hypothesis: True if this is speculative / not yet confirmed.
        direct:     False if interaction is pathway-level, not protein-protein.
    """
    text_refs = {}
    if pmid: text_refs['PMID'] = str(pmid)
    if doi:  text_refs['DOI']  = doi
    return Evidence(
        text=text,
        source_api='manual',
        pmid=pmid,
        text_refs=text_refs or None,
        annotations={
            'date':       str(date.today()),
            'context':    context or 'exploratory',
            'directness': 'direct' if direct else 'indirect',
        },
        epistemics={'hypothesis': hypothesis},
    )


def ag(name, hgnc=None, up=None):
    """Convenience: create an Agent with optional db_refs."""
    db_refs = {}
    if hgnc:
        db_refs['HGNC'] = str(hgnc)
    if up:
        db_refs['UP'] = up
    return Agent(name, db_refs=db_refs or None)


# ── Summary / filtering ─────────────────────────────────────────────────────

def summarize(path=None):
    """Print a readable summary of the store."""
    stmts = load_store(path)
    print(f'{len(stmts)} statements in store\n')
    by_type = {}
    for s in stmts:
        t = type(s).__name__
        by_type.setdefault(t, []).append(s)
    for t, ss in sorted(by_type.items()):
        print(f'  {t}: {len(ss)}')
    print()
    for s in stmts:
        subj = s.agent_list()[0].name if s.agent_list() else '?'
        obj  = s.agent_list()[-1].name if len(s.agent_list()) > 1 else ''
        txt  = s.evidence[0].text[:70] if s.evidence else ''
        print(f'  [{type(s).__name__:20s}]  {subj:12s} → {obj:12s}  "{txt}"')


def stmts_for(gene, path=None):
    """Return all statements involving a given gene."""
    stmts = load_store(path)
    return [
        s for s in stmts
        if any(a.name == gene for a in s.agent_list() if a)
    ]


# ── Context / interaction browsing ───────────────────────────────────────────

def _agents_from_raw(stmt_dict: dict) -> list[str]:
    """Extract agent names from a raw statement JSON dict."""
    names = []
    if 'members' in stmt_dict:
        for m in stmt_dict['members']:
            if m and m.get('name'):
                names.append(m['name'])
    else:
        for key in ('subj', 'enz', 'obj', 'sub'):
            a = stmt_dict.get(key)
            if a and a.get('name'):
                names.append(a['name'])
    return names


def all_contexts(store_path=None) -> dict[str, int]:
    """Return {context_tag: count} from the statement store."""
    p = pathlib.Path(store_path) if store_path else STORE_PATH
    if not p.exists():
        return {}
    with open(p) as f:
        raw = json.load(f)
    contexts: dict[str, int] = {}
    for s in raw:
        for ev in s.get('evidence', []):
            ctx = (ev.get('annotations') or {}).get('context')
            if ctx:
                contexts[ctx] = contexts.get(ctx, 0) + 1
    return dict(sorted(contexts.items()))


def genes_by_context(context: str, store_path=None):
    """Return genes appearing in statements with the given context tag as a GeneList."""
    from geneinfo.genelist import GeneList
    p = pathlib.Path(store_path) if store_path else STORE_PATH
    if not p.exists():
        return GeneList([])
    with open(p) as f:
        raw = json.load(f)
    genes: set[str] = set()
    for s in raw:
        has_ctx = any(
            context.lower() in ((ev.get('annotations') or {}).get('context') or '').lower()
            for ev in s.get('evidence', [])
        )
        if has_ctx:
            for a in _agents_from_raw(s):
                genes.add(a)
    return GeneList(sorted(genes))


def interactors(gene: str, store_path=None) -> list[dict]:
    """
    Return all genes that interact with *gene* in the statement store.

    Returns a list of dicts:
      [{'gene': 'MAPT', 'type': 'Phosphorylation', 'context': '...', 'refs': [...]}, ...]
    """
    p = pathlib.Path(store_path) if store_path else STORE_PATH
    if not p.exists():
        return []
    with open(p) as f:
        raw = json.load(f)
    results = []
    for s in raw:
        agents = _agents_from_raw(s)
        if gene not in agents:
            continue
        partners = [a for a in agents if a != gene]
        ev0 = s.get('evidence', [{}])[0]
        refs = []
        for ev in s.get('evidence', []):
            tr = ev.get('text_refs') or {}
            if tr.get('PMID'):
                refs.append(f"PMID:{tr['PMID']}")
            elif tr.get('DOI'):
                refs.append(f"DOI:{tr['DOI']}")
        for p_name in partners:
            results.append({
                'gene': p_name,
                'type': s.get('type', '?'),
                'context': (ev0.get('annotations') or {}).get('context', ''),
                'refs': refs,
            })
    return results


# ── Query engine ─────────────────────────────────────────────────────────────
# Reuses path-resolution helpers from gene_registry.

from gene_registry import _resolve_path, _split_candidates, _validate_path, _make_path_filter


def query_statements(intersection: bool = True,
                     store_path=None,
                     gene: str = None,
                     context: str = None,
                     **kwargs) -> list[dict]:
    """Query the statement store, returning matching raw statement dicts.

    All keyword arguments are treated as path-based regex filters.
    Underscore-separated names map to nested dict keys, e.g.:

        evidence_text='kinase'         → stmt['evidence'][*]['text']
        evidence_annotations_context='cAMP'
                                       → stmt['evidence'][*]['annotations']['context']
        subj_name='ADRA2C'             → stmt['subj']['name']
        type='Phospho.*'               → stmt['type']

    The special ``gene`` parameter matches against agent names extracted
    from subj/obj/enz/sub/members fields, supporting regex patterns.

    When a path crosses a list (like evidence), the regex is matched
    against every element — a hit on any element counts as a match.

    Parameters
    ----------
    intersection : bool
        True (default) = AND logic; False = OR logic across filters.
    gene : str, optional
        Regex pattern matched against all agent names in a statement.
    context : str, optional
        Alias for ``evidence_annotations_context``. Regex pattern matched
        against the context annotation on each evidence object.
    **kwargs
        Path-based regex filters (see above).

    Returns a list of raw statement dicts (as stored in statements.json).
    """
    import re

    p = pathlib.Path(store_path) if store_path else STORE_PATH
    if not p.exists():
        return []
    with open(p) as f:
        raw = json.load(f)

    # Alias: context → evidence_annotations_context
    if context is not None:
        kwargs.setdefault('evidence_annotations_context', context)

    str_kwargs = {k: v for k, v in kwargs.items() if isinstance(v, str)}
    # Validate paths against a sample of entries
    if raw and str_kwargs:
        sample = raw[:min(5, len(raw))]
        for kwarg_name in str_kwargs:
            _validate_path(kwarg_name, sample)

    filters = [_make_path_filter(k, v) for k, v in str_kwargs.items()]

    # Special gene kwarg — match against agent names
    if gene:
        gene_pat = re.compile(gene, re.IGNORECASE)
        def _gene_filter(s):
            return any(gene_pat.search(a) for a in _agents_from_raw(s))
        filters.append(_gene_filter)

    # Boolean flags
    if kwargs.get('hypothesis_only'):
        filters.append(_make_path_filter(
            'evidence_epistemics_hypothesis', '^True$'))

    if not filters:
        return raw

    combine = all if intersection else any
    return [s for s in raw if combine(f(s) for f in filters)]
