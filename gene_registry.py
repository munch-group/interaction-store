"""
gene_registry.py
----------------
Sidecar gene/agent attribute store for the INDRA interaction notebook.

The registry is a plain JSON dict:  gene_name (str) → attribute dict.
It is intentionally separate from INDRA statements so that gene-level
metadata (chromosome, group membership, analysis provenance, references)
can be queried and updated independently of the interaction graph.
"""

import json
import pathlib
from typing import Any

REGISTRY_PATH = pathlib.Path('gene_registry.json')


# ── I/O ───────────────────────────────────────────────────────────────────────

def load_registry(path=REGISTRY_PATH) -> dict:
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def save_registry(reg: dict, path=REGISTRY_PATH):
    with open(path, 'w') as f:
        json.dump(reg, f, indent=2, sort_keys=True)
    print(f'Saved registry with {len(reg)} genes → {path}')


# ── Mutation helpers ──────────────────────────────────────────────────────────

def add_gene(name: str, attrs: dict, path=REGISTRY_PATH):
    """Add or fully replace a gene entry."""
    reg = load_registry(path)
    reg[name] = attrs
    save_registry(reg, path)


def update_gene(name: str, attrs: dict, path=REGISTRY_PATH):
    """Merge attrs into an existing entry (shallow merge)."""
    reg = load_registry(path)
    existing = reg.get(name, {})
    existing.update(attrs)
    reg[name] = existing
    save_registry(reg, path)


def add_to_group(name: str, group: str, path=REGISTRY_PATH):
    """Append gene to a gene_group without overwriting other attrs."""
    reg = load_registry(path)
    entry = reg.setdefault(name, {})
    groups = entry.setdefault('groups', [])
    if group not in groups:
        groups.append(group)
    save_registry(reg, path)


def add_reference(name: str, pmid: str = None, doi: str = None,
                  note: str = '', path=REGISTRY_PATH):
    """Append a literature reference to a gene entry."""
    reg = load_registry(path)
    entry = reg.setdefault(name, {})
    refs = entry.setdefault('references', [])
    ref = {}
    if pmid: ref['pmid'] = pmid
    if doi:  ref['doi']  = doi
    if note: ref['note'] = note
    refs.append(ref)
    save_registry(reg, path)


# ── Query helpers ─────────────────────────────────────────────────────────────

def get_group(group_name: str, path=REGISTRY_PATH) -> list[str]:
    """Return all gene names belonging to a group."""
    reg = load_registry(path)
    return [
        name for name, attrs in reg.items()
        if group_name in attrs.get('groups', [])
    ]


def get_by_chromosome(chrom: str, path=REGISTRY_PATH) -> list[str]:
    """Return all genes on a given chromosome ('X', 'Y', 'auto', 'mito')."""
    reg = load_registry(path)
    return [
        name for name, attrs in reg.items()
        if attrs.get('chromosome') == chrom
    ]


def get_rescue_candidates(path=REGISTRY_PATH) -> list[str]:
    """Return genes flagged as rescue candidates."""
    reg = load_registry(path)
    return [
        name for name, attrs in reg.items()
        if attrs.get('rescue_logic') in ('rheostat', 'paralog_backup')
    ]


def enrich_graph(G, name_prop, path=REGISTRY_PATH):
    """
    Copy registry attributes onto graph-tool Graph vertex property maps.
    Adds vertex properties: chromosome, rescue_logic.
    Stores groups as a comma-separated string property.

    Parameters
    ----------
    G : graph_tool.Graph
        The interaction graph (vertices must already have *name_prop*).
    name_prop : graph_tool.VertexPropertyMap (type "string")
        Vertex property map containing gene names.

    Returns the graph with new/updated vertex property maps attached as
    G.vp['chromosome'], G.vp['rescue_logic'], G.vp['groups'].
    """
    reg = load_registry(path)

    if 'chromosome' not in G.vp:
        G.vp['chromosome'] = G.new_vertex_property('string')
    if 'rescue_logic' not in G.vp:
        G.vp['rescue_logic'] = G.new_vertex_property('string')
    if 'groups' not in G.vp:
        G.vp['groups'] = G.new_vertex_property('string')

    for v in G.vertices():
        name = name_prop[v]
        info = reg.get(name, {})
        G.vp['chromosome'][v] = info.get('chromosome', 'unknown')
        G.vp['rescue_logic'][v] = info.get('rescue_logic', 'none')
        G.vp['groups'][v] = ','.join(info.get('groups', []))

    return G


def get_gene_info(name: str, path=REGISTRY_PATH) -> dict | None:
    """Return full registry entry for a gene, or None if not found."""
    reg = load_registry(path)
    return reg.get(name)


def get_all_groups(path=REGISTRY_PATH) -> dict[str, list[str]]:
    """Return {group_name: [gene, ...]} for every group in the registry."""
    reg = load_registry(path)
    groups: dict[str, list[str]] = {}
    for name, attrs in reg.items():
        for g in attrs.get('groups', []):
            groups.setdefault(g, []).append(name)
    return {k: sorted(v) for k, v in sorted(groups.items())}


def get_all_contexts(store_path='statements.json') -> dict[str, int]:
    """Return {context_tag: count} from the statement store."""
    import pathlib as _pl
    p = _pl.Path(store_path)
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


def genes_by_context(context: str, store_path='statements.json'):
    """Return genes appearing in statements with the given context tag as a GeneList."""
    from geneinfo.genelist import GeneList
    import pathlib as _pl
    p = _pl.Path(store_path)
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


def genes_by_group(group_name: str, path=REGISTRY_PATH):
    """Return all genes belonging to a group as a GeneList."""
    from geneinfo.genelist import GeneList
    return GeneList(get_group(group_name, path))


def get_interactors(gene: str, store_path='statements.json') -> list[dict]:
    """
    Return all genes that interact with *gene* in the statement store.

    Returns a list of dicts:
      [{'gene': 'MAPT', 'type': 'Phosphorylation', 'context': '...', 'refs': [...]}, ...]
    """
    import pathlib as _pl
    p = _pl.Path(store_path)
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


def _resolve_path(obj, keys: list[str]) -> list:
    """Walk a nested dict/list structure following *keys*.

    At each step: if the current value is a list, descend into every
    element. Returns a flat list of leaf values reached by the path.
    """
    current = [obj]
    for key in keys:
        nxt = []
        for item in current:
            if isinstance(item, dict):
                val = item.get(key)
                if val is not None:
                    nxt.append(val)
            elif isinstance(item, list):
                for elem in item:
                    if isinstance(elem, dict):
                        val = elem.get(key)
                        if val is not None:
                            nxt.append(val)
        current = nxt
    # Flatten any remaining lists so we get individual strings
    flat = []
    for v in current:
        if isinstance(v, list):
            flat.extend(v)
        else:
            flat.append(v)
    return flat


def _split_candidates(name: str) -> list[list[str]]:
    """Generate all ways to split an underscore-separated name into key paths.

    E.g. 'analysis_origin_source' yields:
      [['analysis', 'origin', 'source'],
       ['analysis', 'origin_source'],
       ['analysis_origin', 'source'],
       ['analysis_origin_source']]
    """
    parts = name.split('_')
    if len(parts) == 1:
        return [parts]
    results = []
    # Use bitmask on the (n-1) split points
    for mask in range(1 << (len(parts) - 1)):
        candidate = []
        current = parts[0]
        for i in range(1, len(parts)):
            if mask & (1 << (i - 1)):
                candidate.append(current)
                current = parts[i]
            else:
                current += '_' + parts[i]
        candidate.append(current)
        results.append(candidate)
    return results


def _make_path_filter(kwarg_name: str, pattern: str):
    """Build a filter function from an underscore-separated path and regex.

    Tries all possible splits of the kwarg name into nested dict keys
    to handle keys that themselves contain underscores (e.g.
    'analysis_origin_source' → ['analysis_origin', 'source']).
    A match on any valid split counts.
    """
    import re
    candidates = _split_candidates(kwarg_name)
    pat = re.compile(pattern, re.IGNORECASE)

    def _filter(obj):
        for keys in candidates:
            values = _resolve_path(obj, keys)
            if any(pat.search(str(v)) for v in values):
                return True
        return False
    return _filter


def query_statements(intersection: bool = True,
                     store_path='statements.json',
                     **kwargs) -> list[dict]:
    """Query the statement store, returning matching raw statement dicts.

    All keyword arguments are treated as path-based regex filters.
    Underscore-separated names map to nested dict keys, e.g.:

        evidence_text='kinase'         → stmt['evidence'][*]['text']
        evidence_annotations_context='cAMP'
                                       → stmt['evidence'][*]['annotations']['context']
        subj_name='ADRA2C'             → stmt['subj']['name']
        type='Phospho.*'               → stmt['type']

    When a path crosses a list (like evidence), the regex is matched
    against every element — a hit on any element counts as a match.

    Parameters
    ----------
    intersection : bool
        True (default) = AND logic; False = OR logic across filters.
    **kwargs
        Path-based regex filters (see above).

    Returns a list of raw statement dicts (as stored in statements.json).
    """
    import pathlib as _pl
    p = _pl.Path(store_path)
    if not p.exists():
        return []
    with open(p) as f:
        raw = json.load(f)

    filters = [_make_path_filter(k, v) for k, v in kwargs.items()
               if isinstance(v, str)]
    # Boolean flags
    if kwargs.get('hypothesis_only'):
        filters.append(_make_path_filter(
            'evidence_epistemics_hypothesis', '^True$'))

    if not filters:
        return raw

    combine = all if intersection else any
    return [s for s in raw if combine(f(s) for f in filters)]


def query_genes(intersection: bool = True,
                path=REGISTRY_PATH, **kwargs) -> dict[str, dict]:
    """Query the gene registry, returning matching entries.

    All keyword arguments are treated as path-based regex filters.
    Underscore-separated names map to nested dict keys, e.g.:

        gene='PRK'                       → match gene name
        chromosome='^X$'                 → entry['chromosome']
        analysis_origin_source='IBDmix'  → entry['analysis_origin']['source']
        groups='cAMP'               → entry['groups'][*]
        notes='rescue'                   → entry['notes']
        rescue_logic='rheostat'          → entry['rescue_logic']
        references_pmid='23640'          → entry['references'][*]['pmid']

    When a path crosses a list (like groups or references), the
    regex is matched against every element.

    The special kwarg ``gene`` matches against the gene name (dict key)
    rather than a nested path.

    Parameters
    ----------
    intersection : bool
        True (default) = AND logic; False = OR logic across filters.
    **kwargs
        Path-based regex filters (see above).

    Returns {gene_name: attrs_dict} for matching entries.
    """
    import re
    reg = load_registry(path)
    items = list(reg.items())

    filters = []
    for kwarg_name, pattern in kwargs.items():
        if not isinstance(pattern, str):
            continue
        if kwarg_name == 'gene':
            pat = re.compile(pattern, re.IGNORECASE)
            filters.append(lambda n, a, _p=pat: bool(_p.search(n)))
        else:
            path_filter = _make_path_filter(kwarg_name, pattern)
            filters.append(lambda n, a, _f=path_filter: _f(a))

    if not filters:
        return dict(items)

    combine = all if intersection else any
    return {n: a for n, a in items if combine(f(n, a) for f in filters)}


def _agents_from_raw(stmt_dict: dict) -> list[str]:
    """Extract agent names from a raw statement JSON dict."""
    names = []
    # Complex uses 'members', others use 'subj'/'obj' or 'enz'/'sub'
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


def summarise(path=REGISTRY_PATH):
    reg = load_registry(path)
    print(f'{len(reg)} genes in registry\n')
    by_chrom = {}
    for name, attrs in reg.items():
        c = attrs.get('chromosome', 'unknown')
        by_chrom.setdefault(c, []).append(name)
    for chrom in ('auto', 'X', 'Y', 'mito', 'unknown'):
        genes = by_chrom.get(chrom, [])
        if genes:
            print(f'  {chrom:8s}: {", ".join(sorted(genes))}')
    # Group membership
    all_groups: dict[str, list] = {}
    for name, attrs in reg.items():
        for g in attrs.get('groups', []):
            all_groups.setdefault(g, []).append(name)
    print()
    for group, members in sorted(all_groups.items()):
        print(f'  [{group}]')
        print(f'    {", ".join(sorted(members))}')
