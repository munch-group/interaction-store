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


from pycirclize import Circos
from pycirclize.utils import ColorCycler, load_eukaryote_example_dataset
from collections import defaultdict
import matplotlib
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

chrom_lengths = {'hg19': {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
                            'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
                            'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
                            'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
                            'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
                            'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566},
                    'hg38': {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
                            'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
                            'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
                            'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
                            'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                            'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}}


def register_shifted_cmap(cmap_name, n_lines=24):
    name_shifted = f'shifted_{cmap_name}'
    n = n_lines // 2
    cmap = matplotlib.colormaps[cmap_name]
    colors = cmap(np.linspace(0, 1, n))
    l = list(range(n))
    i = n//2
    colors = [colors[x] for y in zip(l, l[i:] + l[:i]) for x in y][:n_lines]
    _cmap = LinearSegmentedColormap.from_list(f'shifted_{cmap_name}', colors)
    try:
        matplotlib.colormaps.register(cmap=_cmap)
    except ValueError:
        # already registered, ignore
        pass
    return name_shifted

def circos_plot(stmts, reg, cmap=None, scalings={}, assembly='hg38',
        ideogram_base = 97, ideogram_height = 3, figsize = (8, 8)):

    ColorCycler.set_cmap(cmap)
    sector_lengths = chrom_lengths[assembly].copy()
    sector_lengths['chrX'] *= 10
    circos = Circos(sectors=sector_lengths, space=3)
    chr_names = [s.name for s in circos.sectors]
    colors = ColorCycler.get_color_list(len(chr_names))
    chr_name2color = {name: color for name, color in zip(chr_names, colors)}

    label_genes = set([a.name for st in stmts for a in st.agent_list() if a])
    gene_coordinates = {}
    for name, data in reg.items():
        if name in label_genes and 'coordinates' in data:
            gene_coordinates[name] = [
                data['coordinates'][assembly]['chrom'],
                data['coordinates'][assembly]['start'],
                data['coordinates'][assembly]['end']
            ]
    gene_labels = defaultdict(list)
    for name, (chrom, start, end) in gene_coordinates.items():
        gene_labels[chrom].append([(start+end)/2 * scalings.get(chrom, 1), name])

    for sector in circos.sectors:
        sector.text(sector.name, r=105, size=8, color=chr_name2color[sector.name])
        outer_track = sector.add_track((ideogram_base, ideogram_base+ideogram_height))
        outer_track.axis(fc="#eeeeee")

        for pos, label in gene_labels.get(sector.name, []):
            outer_track.annotate(pos, label, label_size=7,
                                min_r=ideogram_base,
                                max_r=ideogram_base+5*ideogram_height,
                                )
    for st in stmts:
        agent_list = st.agent_list()
        assert len(agent_list) == 2
        a, b = agent_list
        if a.name in gene_coordinates and b.name in gene_coordinates:
            _from, _to = gene_coordinates[a.name], gene_coordinates[b.name]
            _from[1] *= scalings.get(_from[0], 1)
            _from[2] *= scalings.get(_from[0], 1)
            _to[1] *= scalings.get(_to[0], 1)
            _to[2] *= scalings.get(_to[0], 1)
            color = chr_name2color[_from[0]]
            circos.link(_from, _to,
                        color=color,
                        lw=0.5,
                        alpha=1,
                        zorder=0)

    fig = circos.plotfig(figsize=figsize)
