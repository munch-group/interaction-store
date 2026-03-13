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
