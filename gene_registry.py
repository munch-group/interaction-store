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
    groups = entry.setdefault('gene_groups', [])
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
        if group_name in attrs.get('gene_groups', [])
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
    Stores gene_groups as a comma-separated string property.

    Parameters
    ----------
    G : graph_tool.Graph
        The interaction graph (vertices must already have *name_prop*).
    name_prop : graph_tool.VertexPropertyMap (type "string")
        Vertex property map containing gene names.

    Returns the graph with new/updated vertex property maps attached as
    G.vp['chromosome'], G.vp['rescue_logic'], G.vp['gene_groups'].
    """
    reg = load_registry(path)

    if 'chromosome' not in G.vp:
        G.vp['chromosome'] = G.new_vertex_property('string')
    if 'rescue_logic' not in G.vp:
        G.vp['rescue_logic'] = G.new_vertex_property('string')
    if 'gene_groups' not in G.vp:
        G.vp['gene_groups'] = G.new_vertex_property('string')

    for v in G.vertices():
        name = name_prop[v]
        info = reg.get(name, {})
        G.vp['chromosome'][v] = info.get('chromosome', 'unknown')
        G.vp['rescue_logic'][v] = info.get('rescue_logic', 'none')
        G.vp['gene_groups'][v] = ','.join(info.get('gene_groups', []))

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
        for g in attrs.get('gene_groups', []):
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


def query_statements(gene: str = None, stmt_type: str = None,
                     context: str = None, hypothesis_only: bool = False,
                     store_path='statements.json') -> list[dict]:
    """Query the statement store, returning matching raw statement dicts.

    Parameters
    ----------
    gene : str, optional
        Filter to statements involving this gene (subject or object).
    stmt_type : str, optional
        Filter by INDRA statement type name (e.g. 'Phosphorylation').
    context : str, optional
        Case-insensitive substring match on evidence context annotation.
    hypothesis_only : bool
        If True, return only statements with hypothesis evidence.

    Returns a list of raw statement dicts (as stored in statements.json).
    """
    import pathlib as _pl
    p = _pl.Path(store_path)
    if not p.exists():
        return []
    with open(p) as f:
        raw = json.load(f)
    results = raw
    if gene:
        results = [s for s in results if gene in _agents_from_raw(s)]
    if stmt_type:
        results = [s for s in results if s.get('type') == stmt_type]
    if context:
        results = [s for s in results
                   if any(context.lower() in
                          ((ev.get('annotations') or {}).get('context') or '').lower()
                          for ev in s.get('evidence', []))]
    if hypothesis_only:
        results = [s for s in results
                   if any((ev.get('epistemics') or {}).get('hypothesis')
                          for ev in s.get('evidence', []))]
    return results


def query_registry(gene: str = None, group: str = None,
                   chromosome: str = None, rescue_only: bool = False,
                   analysis: str = None,
                   path=REGISTRY_PATH) -> dict[str, dict]:
    """Query the gene registry, returning matching entries.

    Parameters
    ----------
    gene : str, optional
        Return only this gene's entry.
    group : str, optional
        Filter to genes in this group.
    chromosome : str, optional
        Filter by chromosome ('auto', 'X', 'Y', 'mito').
    rescue_only : bool
        If True, return only genes with non-trivial rescue_logic.
    analysis : str, optional
        Substring match on analysis_origin.analysis.

    Returns {gene_name: attrs_dict} for matching entries.
    """
    reg = load_registry(path)
    items = list(reg.items())
    if gene:
        items = [(n, a) for n, a in items if n == gene]
    if group:
        items = [(n, a) for n, a in items
                 if group in a.get('gene_groups', [])]
    if chromosome:
        items = [(n, a) for n, a in items
                 if a.get('chromosome') == chromosome]
    if rescue_only:
        items = [(n, a) for n, a in items
                 if a.get('rescue_logic') not in (None, 'none')]
    if analysis:
        items = [(n, a) for n, a in items
                 if analysis in (a.get('analysis_origin') or {}).get('analysis', '')]
    return dict(items)


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
        for g in attrs.get('gene_groups', []):
            all_groups.setdefault(g, []).append(name)
    print()
    for group, members in sorted(all_groups.items()):
        print(f'  [{group}]')
        print(f'    {", ".join(sorted(members))}')
