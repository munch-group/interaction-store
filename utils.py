"""
utils.py
--------
Notebook-specific helpers: graph-tool traversal and circos visualization.

Statement I/O and query functions live in statement_store.py.
"""

import graph_tool.all as gt


def _name_to_vertex(G):
    """Build a name→vertex lookup from G.vp['name']."""
    return {G.vp['name'][v]: v for v in G.vertices()}

def neighbors(gene, G=None):
    """Print all direct interaction partners of a gene."""
    if G is None:
        G = globals()['G']
    lookup = _name_to_vertex(G)
    if gene not in lookup:
        print(f'{gene} not in graph')
        return
    v = lookup[gene]

    print(f'\n{gene} → (outgoing)')
    for e in v.out_edges():
        tgt = G.vp['name'][e.target()]
        st  = G.ep['stmt_type'][e]
        print(f'  {st:20s} → {tgt}')

    print(f'\n{gene} ← (incoming)')
    for e in v.in_edges():
        src = G.vp['name'][e.source()]
        st  = G.ep['stmt_type'][e]
        print(f'  {src:20s} → [{st}]')


def path_between(source, target, G=None):
    """Find shortest path between two genes."""
    if G is None:
        G = globals()['G']
    lookup = _name_to_vertex(G)
    if source not in lookup:
        print(f'Node not found: {source}')
        return
    if target not in lookup:
        print(f'Node not found: {target}')
        return
    vs = lookup[source]
    vt = lookup[target]
    vlist, elist = gt.shortest_path(G, vs, vt)
    if not vlist:
        print(f'No path from {source} to {target}')
    else:
        names = [G.vp['name'][v] for v in vlist]
        print(' → '.join(names))

