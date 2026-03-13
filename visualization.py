"""
visualization.py
----------------
Interactive (ipycytoscape) and static (graph-tool) network visualization
for the INDRA interaction store.

Public API
----------
interactive_network(stmts, reg, G=None, highlight=None) -> ipywidgets.VBox
save_static_png(stmts, reg, G=None, output='interaction_network.png', size=(1800,1300)) -> str
"""

# ── colour / style constants ────────────────────────────────────────────────

MODULE_PRIORITY = [
    'cAMP/PKA module',
    'IRS2/Akt module',
    'MT lattice/transport',
    'Neurotrophin/endosome',
    'NMD/RNA processing',
    'NF-kB signalling',
    'Centriolar/manchette',
    'Phosphoinositide',
]

MODULE_BORDER_COLORS = {
    'cAMP/PKA module':       '#8E44AD',
    'IRS2/Akt module':       '#2980B9',
    'MT lattice/transport':  '#27AE60',
    'Neurotrophin/endosome': '#E67E22',
    'NMD/RNA processing':    '#1ABC9C',
    'NF-kB signalling':      '#E74C3C',
    'Centriolar/manchette':  '#F39C12',
    'Phosphoinositide':      '#3498DB',
}

CHROM_BG_COLORS = {
    'auto':    '#5B9BD5',
    'X':       '#D95F5F',
    'Y':       '#E8A020',
    'unknown': '#DADDE1',
}

RGBA_NODE = {
    'auto':    [0.357, 0.608, 0.835, 1.0],
    'X':       [0.851, 0.373, 0.373, 1.0],
    'Y':       [0.910, 0.627, 0.125, 1.0],
    'unknown': [0.855, 0.867, 0.882, 1.0],
}

RGBA_EDGE = {
    'Activation':      [0.153, 0.682, 0.376, 0.8],
    'Inhibition':      [0.753, 0.224, 0.169, 0.8],
    'Phosphorylation': [0.161, 0.502, 0.725, 0.8],
    'IncreaseAmount':  [0.153, 0.682, 0.376, 0.8],
    'DecreaseAmount':  [0.753, 0.224, 0.169, 0.8],
    'Complex':         [0.533, 0.533, 0.533, 0.8],
}
DEFAULT_RGBA_EDGE = [0.667, 0.667, 0.667, 0.8]


# ── private helpers ──────────────────────────────────────────────────────────

def _primary_module(gene, reg):
    """First matching module from registry gene_groups."""
    gene_groups = reg.get(gene, {}).get('gene_groups', [])
    for mod in MODULE_PRIORITY:
        if mod in gene_groups:
            return mod
    return 'none'


def _build_edge_evidence(stmts):
    """Build (src, tgt, type) -> (ev_text, ctx) dict from statements."""
    edge_evidence = {}
    for s in stmts:
        agents = [a for a in s.agent_list() if a is not None]
        stype = type(s).__name__
        if stype == 'Complex' and len(agents) >= 2:
            for i, a in enumerate(agents):
                for b in agents[i + 1:]:
                    key = (a.name, b.name, 'Complex')
                    ev_text = s.evidence[0].text[:200] if s.evidence else ''
                    ctx = s.evidence[0].annotations.get('context', '') if s.evidence else ''
                    edge_evidence[key] = (ev_text, ctx)
        elif len(agents) >= 2:
            key = (agents[0].name, agents[-1].name, stype)
            ev_text = s.evidence[0].text[:200] if s.evidence else ''
            ctx = s.evidence[0].annotations.get('context', '') if s.evidence else ''
            edge_evidence[key] = (ev_text, ctx)
    return edge_evidence


def _build_graph(stmts, reg):
    """Build a graph-tool directed Graph from INDRA statements + registry."""
    import graph_tool.all as gt

    G = gt.Graph(directed=True)
    G.vp['name'] = G.new_vertex_property('string')
    G.vp['chromosome'] = G.new_vertex_property('string')
    G.ep['stmt_type'] = G.new_edge_property('string')

    node_idx = {}

    def _get_or_add(name):
        if name not in node_idx:
            v = G.add_vertex()
            node_idx[name] = int(v)
            G.vp['name'][v] = name
            info = reg.get(name, {})
            G.vp['chromosome'][v] = info.get('chromosome', 'unknown')
        return G.vertex(node_idx[name])

    for s in stmts:
        agents = [a for a in s.agent_list() if a is not None]
        stype = type(s).__name__
        if stype == 'Complex':
            for i, a in enumerate(agents):
                for b in agents[i + 1:]:
                    va, vb = _get_or_add(a.name), _get_or_add(b.name)
                    e = G.add_edge(va, vb)
                    G.ep['stmt_type'][e] = 'Complex'
        elif len(agents) >= 2:
            va = _get_or_add(agents[0].name)
            vb = _get_or_add(agents[-1].name)
            e = G.add_edge(va, vb)
            G.ep['stmt_type'][e] = stype

    return G


def _cytoscape_json(G, reg, edge_evidence):
    """Convert graph-tool Graph to {nodes, edges} dict for ipycytoscape."""
    nodes = []
    edges = []

    for v in G.vertices():
        name = G.vp['name'][v]
        chrom = G.vp['chromosome'][v]
        mod = _primary_module(name, reg)
        rescue = reg.get(name, {}).get('rescue_logic', 'none')
        groups = ', '.join(reg.get(name, {}).get('gene_groups', []))
        info = f'{name}\nChromosome: {chrom}\nModule: {mod}\nGroups: {groups}\nRescue: {rescue}'
        nodes.append({
            'data': {
                'id': name,
                'label': name,
                'chromosome': chrom,
                'module': mod,
                '_info': info,
            },
        })

    for e in G.edges():
        src = G.vp['name'][e.source()]
        tgt = G.vp['name'][e.target()]
        st = G.ep['stmt_type'][e]
        ev_text, ctx = edge_evidence.get((src, tgt, st), ('', ''))
        info = f'{st}: {src} → {tgt}\nContext: {ctx}\nEvidence: {ev_text}'
        edges.append({
            'data': {
                'source': src,
                'target': tgt,
                'stmt_type': st,
                '_info': info,
            },
            'classes': st + ' directed',
        })

    return {'nodes': nodes, 'edges': edges}


def _default_style():
    """Return the CSS selector list for ipycytoscape."""
    style = [
        # Default node
        {'selector': 'node',
         'style': {
             'label': 'data(label)',
             'text-valign': 'center',
             'text-halign': 'center',
             'font-size': '10px',
             'font-weight': 'bold',
             'color': 'white',
             'text-outline-color': '#333',
             'text-outline-width': 1,
             'width': 45,
             'height': 45,
             'background-color': '#DADDE1',
             'border-width': 2,
             'border-color': '#999',
         }},
    ]

    # Chromosome colouring
    for chrom, color in CHROM_BG_COLORS.items():
        if chrom == 'unknown':
            continue
        style.append({
            'selector': f'node[chromosome = "{chrom}"]',
            'style': {'background-color': color},
        })

    # Module border colours
    for mod, color in MODULE_BORDER_COLORS.items():
        style.append({
            'selector': f'node[module = "{mod}"]',
            'style': {'border-color': color},
        })

    style += [
        # Default edge
        {'selector': 'edge',
         'style': {
             'curve-style': 'bezier',
             'target-arrow-shape': 'triangle',
             'target-arrow-color': '#AAAAAA',
             'line-color': '#AAAAAA',
             'width': 2,
             'arrow-scale': 1.2,
         }},
        # Activation — green arrow
        {'selector': '.Activation, .IncreaseAmount',
         'style': {
             'line-color': '#27AE60',
             'target-arrow-color': '#27AE60',
             'target-arrow-shape': 'triangle',
         }},
        # Inhibition — red flat-head
        {'selector': '.Inhibition, .DecreaseAmount',
         'style': {
             'line-color': '#C0392B',
             'target-arrow-color': '#C0392B',
             'target-arrow-shape': 'tee',
         }},
        # Phosphorylation — blue diamond
        {'selector': '.Phosphorylation',
         'style': {
             'line-color': '#2980B9',
             'target-arrow-color': '#2980B9',
             'target-arrow-shape': 'diamond',
         }},
        # Complex — grey dashed, no arrow
        {'selector': '.Complex',
         'style': {
             'line-color': '#888888',
             'line-style': 'dashed',
             'target-arrow-shape': 'none',
         }},
        # Selection highlight
        {'selector': 'node:selected',
         'style': {
             'border-width': 4,
             'border-color': '#FFD700',
         }},
        {'selector': 'edge:selected',
         'style': {
             'width': 4,
             'line-color': '#FFD700',
             'target-arrow-color': '#FFD700',
         }},
    ]
    return style


def _make_click_handlers(info_box):
    """Return (on_node_click, on_edge_click) closures."""
    _style = (
        'padding:8px; font-size:13px; color:#333; '
        'background:#f8f8f8; border:1px solid #ddd; border-radius:4px; '
        'font-family:monospace; white-space:pre-wrap;'
    )

    def on_node_click(event):
        node = event['data']
        tip = node.get('_info', node.get('id', ''))
        html = tip.replace('\n', '<br>')
        info_box.value = (
            f'<div style="{_style}">'
            f'<b>Node:</b> {node.get("id","")}<br>{html}</div>'
        )

    def on_edge_click(event):
        edge = event['data']
        tip = edge.get('_info', '')
        html = tip.replace('\n', '<br>')
        info_box.value = (
            f'<div style="{_style}">'
            f'<b>Edge:</b> {edge.get("source","")} → {edge.get("target","")}<br>{html}</div>'
        )

    return on_node_click, on_edge_click


# ── public API ───────────────────────────────────────────────────────────────

def interactive_network(stmts, reg, G=None, highlight=None):
    """
    Build an interactive ipycytoscape widget with click-to-inspect info box.

    Parameters
    ----------
    stmts : list[Statement]
        INDRA statements (used for edge evidence text).
    reg : dict
        Gene registry (gene_name -> attributes).
    G : graph_tool.Graph, optional
        Pre-built graph. If None, built from stmts + reg.
    highlight : str, optional
        Gene name to highlight (not yet implemented in interactive view).

    Returns
    -------
    ipywidgets.VBox containing the cytoscape widget and info box.
    """
    import ipycytoscape
    import ipywidgets as widgets

    if G is None:
        G = _build_graph(stmts, reg)

    edge_evidence = _build_edge_evidence(stmts)
    cyto_json = _cytoscape_json(G, reg, edge_evidence)

    info_box = widgets.HTML(
        value='<div style="padding:8px; font-size:13px; color:#333; '
              'background:#f8f8f8; border:1px solid #ddd; border-radius:4px; '
              'min-height:40px; font-family:monospace; white-space:pre-wrap;">'
              'Click a node or edge to see details.</div>',
    )

    cyto = ipycytoscape.CytoscapeWidget()
    cyto.graph.add_graph_from_json(cyto_json, directed=True)
    cyto.set_style(_default_style())
    cyto.set_layout(
        name='cose', nodeOverlap=20, idealEdgeLength=80,
        edgeElasticity=100, nestingFactor=1.2,
        gravity=0.5, numIter=1000,
        initialTemp=200, coolingFactor=0.95,
        minTemp=1.0, animate=False,
    )

    on_node_click, on_edge_click = _make_click_handlers(info_box)
    cyto.on('node', 'click', on_node_click)
    cyto.on('edge', 'click', on_edge_click)

    return widgets.VBox([cyto, info_box])


def save_static_png(stmts, reg, G=None, output='interaction_network.png',
                    size=(1800, 1300)):
    """
    Render a static PNG of the interaction network via graph-tool.

    Parameters
    ----------
    stmts : list[Statement]
        INDRA statements.
    reg : dict
        Gene registry.
    G : graph_tool.Graph, optional
        Pre-built graph. If None, built from stmts + reg.
    output : str
        Output file path.
    size : tuple[int, int]
        Image dimensions in pixels.

    Returns
    -------
    str — output file path.
    """
    import graph_tool.all as gt

    if G is None:
        G = _build_graph(stmts, reg)

    # Vertex fill colour from chromosome
    vfill = G.new_vertex_property('vector<double>')
    for v in G.vertices():
        chrom = G.vp['chromosome'][v]
        vfill[v] = RGBA_NODE.get(chrom, RGBA_NODE['unknown'])

    # Edge colour and dash
    ecolor = G.new_edge_property('vector<double>')
    edash = G.new_edge_property('vector<double>')
    for e in G.edges():
        st = G.ep['stmt_type'][e]
        ecolor[e] = RGBA_EDGE.get(st, DEFAULT_RGBA_EDGE)
        edash[e] = [0.02, 0.01] if st == 'Complex' else []

    pos = gt.sfdp_layout(G)
    gt.graph_draw(
        G, pos=pos,
        vertex_fill_color=vfill,
        vertex_color=[0.2, 0.2, 0.2, 0.6],
        vertex_size=28,
        vertex_text=G.vp['name'],
        vertex_text_color=[1, 1, 1, 1],
        vertex_font_size=8,
        vertex_font_weight=1,
        edge_color=ecolor,
        edge_pen_width=2.0,
        edge_dash_style=edash,
        edge_marker_size=12,
        output_size=size,
        output=output,
    )
    print(f'Static PNG → {output}')
    return output
