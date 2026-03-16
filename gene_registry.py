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

REGISTRY_PATH = pathlib.Path(__file__).parent / 'agents.json'


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

def load_registry(path=REGISTRY_PATH) -> dict:
    if isinstance(path, str):
        path = pathlib.Path(path)
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def save_registry(reg: dict, path=REGISTRY_PATH):
    """Save registry atomically."""
    path = pathlib.Path(path)
    _atomic_json_write(path, reg, sort_keys=True)
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

def group_members(group_name: str, path=REGISTRY_PATH) -> list[str]:
    """Return all gene names belonging to a group."""
    reg = load_registry(path)
    return [
        name for name, attrs in reg.items()
        if group_name in attrs.get('groups', [])
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


def gene_data(name: str, path=REGISTRY_PATH) -> dict | None:
    """Return full registry entry for a gene, or None if not found."""
    reg = load_registry(path)
    return reg.get(name)


def all_groups(path=REGISTRY_PATH) -> dict[str, list[str]]:
    """Return {group_name: [gene, ...]} for every group in the registry."""
    reg = load_registry(path)
    groups: dict[str, list[str]] = {}
    for name, attrs in reg.items():
        for g in attrs.get('groups', []):
            groups.setdefault(g, []).append(name)
    return {k: sorted(v) for k, v in sorted(groups.items())}


def genes_by_group(group_name: str, path=REGISTRY_PATH):
    """Return all genes belonging to a group as a GeneList."""
    from geneinfo.genelist import GeneList
    return GeneList(group_members(group_name, path))


# ── Query engine ─────────────────────────────────────────────────────────────

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


def _validate_path(kwarg_name: str, sample_objects: list[dict]):
    """Check that at least one split of kwarg_name resolves to data in
    at least one sample object. Raises KeyError if no split works."""
    candidates = _split_candidates(kwarg_name)
    for obj in sample_objects:
        for keys in candidates:
            if _resolve_path(obj, keys):
                return  # valid
    valid_splits = [' → '.join(c) for c in candidates]
    raise KeyError(
        f"'{kwarg_name}' does not match any field. "
        f"Tried paths: {', '.join(valid_splits)}"
    )


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

    str_kwargs = {k: v for k, v in kwargs.items() if isinstance(v, str)}
    # Validate paths against a sample of entries (skip 'gene' — it matches keys)
    if items and str_kwargs:
        sample = [a for _, a in items[:min(5, len(items))]]
        for kwarg_name in str_kwargs:
            if kwarg_name != 'gene':
                _validate_path(kwarg_name, sample)

    filters = []
    for kwarg_name, pattern in str_kwargs.items():
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
