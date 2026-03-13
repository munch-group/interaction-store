"""
visualization.py
----------------
Interactive (ipycytoscape) and static (graph-tool) network visualization
for the INDRA interaction store.

Public API
----------
interactive_network(stmts, reg, G=None, highlight=None, **kwargs) -> ipywidgets.VBox
save_static_png(stmts, reg, G=None, output='interaction_network.png', size=(1800,1300), **kwargs) -> str

Style kwargs (all optional, applied to both functions where applicable)
-----------------------------------------------------------------------
chrom_colors     : dict  — chromosome -> hex colour, e.g. {'auto': '#5B9BD5', 'X': '#D95F5F'}
module_colors    : dict  — module name -> hex border colour
edge_colors      : dict  — statement type -> hex colour, e.g. {'Activation': '#27AE60'}
node_size        : int   — node diameter in px (default 45 interactive, 28 static)
node_shape       : str   — cytoscape.js shape name (default 'ellipse', which renders as a
                   circle when width==height). Other options: 'rectangle', 'roundrectangle',
                   'diamond', 'triangle', 'pentagon', 'hexagon', 'heptagon', 'octagon',
                   'star', 'barrel', 'tag', 'rhomboid', 'vee'. Static PNG always circles.
font_size        : int|str — label font size (default '10px' interactive, 8 static)
text_color       : str   — node label colour (default 'white')
font_weight      : str   — node label weight: 'bold', 'normal', etc. (default 'bold')
text_outline     : bool  — dark outline around label text (default True)
border_width     : int   — node border width in px (default 2)
edge_width       : int   — edge line width in px (default 2)
background_color : str   — canvas/background colour (default 'transparent');
                   for static PNG default is '#FFFFFF'
"""

# ── default colour / style constants ─────────────────────────────────────────

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

DEFAULT_MODULE_BORDER_COLORS = {
    'cAMP/PKA module':       '#8E44AD',
    'IRS2/Akt module':       '#2980B9',
    'MT lattice/transport':  '#27AE60',
    'Neurotrophin/endosome': '#E67E22',
    'NMD/RNA processing':    '#1ABC9C',
    'NF-kB signalling':      '#E74C3C',
    'Centriolar/manchette':  '#F39C12',
    'Phosphoinositide':      '#3498DB',
}

DEFAULT_CHROM_BG_COLORS = {
    'auto':    '#5B9BD5',
    'X':       '#D95F5F',
    'Y':       '#E8A020',
    'unknown': '#DADDE1',
}

DEFAULT_EDGE_COLORS = {
    'Activation':      '#27AE60',
    'Inhibition':      '#C0392B',
    'Phosphorylation': '#2980B9',
    'IncreaseAmount':  '#27AE60',
    'DecreaseAmount':  '#C0392B',
    'Complex':         '#888888',
}

# RGBA versions for graph-tool static PNG
DEFAULT_RGBA_NODE = {
    'auto':    [0.357, 0.608, 0.835, 1.0],
    'X':       [0.851, 0.373, 0.373, 1.0],
    'Y':       [0.910, 0.627, 0.125, 1.0],
    'unknown': [0.855, 0.867, 0.882, 1.0],
}

DEFAULT_RGBA_EDGE = {
    'Activation':      [0.153, 0.682, 0.376, 0.8],
    'Inhibition':      [0.753, 0.224, 0.169, 0.8],
    'Phosphorylation': [0.161, 0.502, 0.725, 0.8],
    'IncreaseAmount':  [0.153, 0.682, 0.376, 0.8],
    'DecreaseAmount':  [0.753, 0.224, 0.169, 0.8],
    'Complex':         [0.533, 0.533, 0.533, 0.8],
}
_FALLBACK_RGBA_EDGE = [0.667, 0.667, 0.667, 0.8]

EDGE_ARROW_SHAPES = {
    'Activation':      'triangle',
    'IncreaseAmount':  'triangle',
    'Inhibition':      'tee',
    'DecreaseAmount':  'tee',
    'Phosphorylation': 'diamond',
    'Complex':         'none',
}


# ── private helpers ──────────────────────────────────────────────────────────

def _hex_to_rgba(hex_color, alpha=1.0):
    """Convert '#RRGGBB' to [r, g, b, a] in 0-1 range."""
    h = hex_color.lstrip('#')
    return [int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4)] + [alpha]


def _primary_module(gene, reg):
    """First matching module from registry groups."""
    groups = reg.get(gene, {}).get('groups', [])
    for mod in MODULE_PRIORITY:
        if mod in groups:
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
        groups = ', '.join(reg.get(name, {}).get('groups', []))
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
        info = f'{st}: {src} \u2192 {tgt}\nContext: {ctx}\nEvidence: {ev_text}'
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


def _build_style(chrom_colors, module_colors, edge_colors,
                 node_size, node_shape, font_size, border_width, edge_width,
                 text_color, font_weight, text_outline):
    """Return the CSS selector list for ipycytoscape."""
    node_style = {
        'label': 'data(label)',
        'text-valign': 'center',
        'text-halign': 'center',
        'font-size': f'{font_size}px' if isinstance(font_size, int) else font_size,
        'font-weight': font_weight,
        'color': text_color,
        'width': node_size,
        'height': node_size,
        'shape': node_shape,
        'background-color': chrom_colors.get('unknown', '#DADDE1'),
        'border-width': border_width,
        'border-color': '#999',
    }
    if text_outline:
        node_style['text-outline-color'] = '#333'
        node_style['text-outline-width'] = 1
    else:
        node_style['text-outline-width'] = 0

    style = [
        {'selector': 'node', 'style': node_style},
    ]

    # Chromosome colouring
    for chrom, color in chrom_colors.items():
        if chrom == 'unknown':
            continue
        style.append({
            'selector': f'node[chromosome = "{chrom}"]',
            'style': {'background-color': color},
        })

    # Module border colours
    for mod, color in module_colors.items():
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
             'width': edge_width,
             'arrow-scale': 1.2,
         }},
        # Activation — green arrow
        {'selector': '.Activation, .IncreaseAmount',
         'style': {
             'line-color': edge_colors.get('Activation', '#27AE60'),
             'target-arrow-color': edge_colors.get('Activation', '#27AE60'),
             'target-arrow-shape': 'triangle',
         }},
        # Inhibition — red flat-head
        {'selector': '.Inhibition, .DecreaseAmount',
         'style': {
             'line-color': edge_colors.get('Inhibition', '#C0392B'),
             'target-arrow-color': edge_colors.get('Inhibition', '#C0392B'),
             'target-arrow-shape': 'tee',
         }},
        # Phosphorylation — blue diamond
        {'selector': '.Phosphorylation',
         'style': {
             'line-color': edge_colors.get('Phosphorylation', '#2980B9'),
             'target-arrow-color': edge_colors.get('Phosphorylation', '#2980B9'),
             'target-arrow-shape': 'diamond',
         }},
        # Complex — grey dashed, no arrow
        {'selector': '.Complex',
         'style': {
             'line-color': edge_colors.get('Complex', '#888888'),
             'line-style': 'dashed',
             'target-arrow-shape': 'none',
         }},
        # Selection highlight
        {'selector': 'node:selected',
         'style': {
             'border-width': border_width + 2,
             'border-color': '#FFD700',
         }},
        {'selector': 'edge:selected',
         'style': {
             'width': edge_width + 2,
             'line-color': '#FFD700',
             'target-arrow-color': '#FFD700',
         }},
    ]
    return style


def _build_legend_html(chrom_colors, module_colors, edge_colors):
    """Build an HTML legend string for the info box."""
    swatch = (
        '<span style="display:inline-block; width:14px; height:14px; '
        'border-radius:50%; background:{color}; border:2px solid {border}; '
        'vertical-align:middle; margin-right:6px;"></span>'
    )
    line_swatch = (
        '<span style="display:inline-block; width:24px; height:0; '
        'border-top:3px {dash} {color}; vertical-align:middle; '
        'margin-right:6px;"></span>'
    )
    arrow_label = {
        'Activation':      '\u25b6 triangle',
        'IncreaseAmount':  '\u25b6 triangle',
        'Inhibition':      '\u22a3 flat-head',
        'DecreaseAmount':  '\u22a3 flat-head',
        'Phosphorylation': '\u25c6 diamond',
        'Complex':         'none (dashed)',
    }

    parts = ['<b>Legend</b><br><br>']

    # Node fill = chromosome
    parts.append('<b>Node fill \u2014 Chromosome</b><br>')
    chrom_labels = {'auto': 'Autosomal', 'X': 'X-linked', 'Y': 'Y-linked', 'unknown': 'Unknown/intermediate'}
    for chrom, color in chrom_colors.items():
        label = chrom_labels.get(chrom, chrom)
        parts.append(swatch.format(color=color, border='#666') + f'{label}<br>')

    # Node border = module
    parts.append('<br><b>Node border \u2014 Primary module</b><br>')
    for mod, color in module_colors.items():
        parts.append(swatch.format(color='#ccc', border=color) + f'{mod}<br>')

    # Edges
    parts.append('<br><b>Edges \u2014 Interaction type</b><br>')
    seen = set()
    for stype in ['Activation', 'Inhibition', 'Phosphorylation',
                   'IncreaseAmount', 'DecreaseAmount', 'Complex']:
        color = edge_colors.get(stype, '#AAAAAA')
        dash = 'dashed' if stype == 'Complex' else 'solid'
        arrow = arrow_label.get(stype, 'triangle')
        label = f'{stype} \u2014 {arrow}'
        # Merge Activation/IncreaseAmount and Inhibition/DecreaseAmount
        key = (color, arrow)
        if key in seen:
            continue
        seen.add(key)
        if stype == 'IncreaseAmount':
            label = 'Activation / IncreaseAmount \u2014 ' + arrow
        elif stype == 'DecreaseAmount':
            label = 'Inhibition / DecreaseAmount \u2014 ' + arrow
        parts.append(line_swatch.format(color=color, dash=dash) + f'{label}<br>')

    parts.append('<br><i>Click a node or edge for details.</i>')
    return ''.join(parts)


def _make_click_handlers(info_box, legend_html):
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
            f'<b>Edge:</b> {edge.get("source","")} \u2192 {edge.get("target","")}<br>{html}</div>'
        )

    def on_legend_click(_):
        info_box.value = (
            f'<div style="{_style}">{legend_html}</div>'
        )

    return on_node_click, on_edge_click, on_legend_click


# ── public API ───────────────────────────────────────────────────────────────

def interactive_network(stmts, reg, G=None, highlight=None, **kwargs):
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

    Style kwargs (all optional)
    ---------------------------
    chrom_colors     : dict — chromosome -> hex colour
    module_colors    : dict — module name -> hex border colour
    edge_colors      : dict — statement type -> hex colour
    node_size        : int  — node diameter in px (default 45)
    node_shape       : str  — cytoscape.js shape: 'ellipse' (circle, default),
                       'rectangle', 'roundrectangle', 'diamond', 'triangle',
                       'pentagon', 'hexagon', 'star', 'tag', 'vee', etc.
    font_size        : int|str — label font size (default '10px')
    text_color       : str  — node label colour (default 'white')
    font_weight      : str  — label weight: 'bold', 'normal' (default 'bold')
    text_outline     : bool — dark outline around labels (default True)
    border_width     : int  — node border width (default 2)
    edge_width       : int  — edge line width (default 2)
    background_color : str  — canvas background (default 'transparent')

    Returns
    -------
    ipywidgets.VBox containing the cytoscape widget, legend button, and info box.
    """
    import ipycytoscape
    import ipywidgets as widgets
    from IPython.display import display, HTML as DisplayHTML

    if G is None:
        G = _build_graph(stmts, reg)

    # Resolve style kwargs with defaults
    chrom_colors = {**DEFAULT_CHROM_BG_COLORS, **(kwargs.get('chrom_colors') or {})}
    module_colors = {**DEFAULT_MODULE_BORDER_COLORS, **(kwargs.get('module_colors') or {})}
    edge_colors = {**DEFAULT_EDGE_COLORS, **(kwargs.get('edge_colors') or {})}
    node_size = kwargs.get('node_size', 45)
    node_shape = kwargs.get('node_shape', 'ellipse')
    font_size = kwargs.get('font_size', '10px')
    border_width = kwargs.get('border_width', 2)
    edge_width = kwargs.get('edge_width', 2)
    background_color = kwargs.get('background_color', 'transparent')
    text_color = kwargs.get('text_color', 'white')
    font_weight = kwargs.get('font_weight', 'bold')
    text_outline = kwargs.get('text_outline', True)

    edge_evidence = _build_edge_evidence(stmts)
    cyto_json = _cytoscape_json(G, reg, edge_evidence)

    _box_style = (
        'padding:8px; font-size:13px; color:#333; '
        'background:#f8f8f8; border:1px solid #ddd; border-radius:4px; '
        'min-height:40px; font-family:monospace; white-space:pre-wrap;'
    )
    legend_html = _build_legend_html(chrom_colors, module_colors, edge_colors)
    info_box = widgets.HTML(
        value=f'<div style="{_box_style}">Click a node or edge for details, or press Show Legend.</div>',
    )

    cyto = ipycytoscape.CytoscapeWidget()
    cyto.graph.add_graph_from_json(cyto_json, directed=True)
    cyto.set_style(_build_style(
        chrom_colors, module_colors, edge_colors,
        node_size, node_shape, font_size, border_width, edge_width,
        text_color, font_weight, text_outline,
    ))
    cyto.set_layout(
        name='cose', nodeOverlap=20, idealEdgeLength=80,
        edgeElasticity=100, nestingFactor=1.2,
        gravity=0.5, numIter=1000,
        initialTemp=200, coolingFactor=0.95,
        minTemp=1.0, animate=False,
    )
    cyto.layout.height = '600px'

    # Override ipycytoscape's .custom-widget background via IPython.display.
    # This injects CSS into the cell output context before the widget renders.
    from IPython.display import display as ipy_display, HTML as IPyHTML
    if background_color != 'transparent':
        ipy_display(IPyHTML(
            '<style>'
            f'.custom-widget {{ background-color: {background_color} !important; }}'
            '</style>'
        ))

    on_node_click, on_edge_click, on_legend_click = _make_click_handlers(info_box, legend_html)
    cyto.on('node', 'click', on_node_click)
    cyto.on('edge', 'click', on_edge_click)

    legend_btn = widgets.Button(
        description='Show Legend',
        button_style='info',
        layout=widgets.Layout(width='120px', margin='4px 0'),
    )
    legend_btn.on_click(on_legend_click)

    return widgets.VBox([cyto, legend_btn, info_box])


def save_static_png(stmts, reg, G=None, output='interaction_network.png',
                    size=(1800, 1300), **kwargs):
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

    Style kwargs (all optional)
    ---------------------------
    chrom_colors     : dict — chromosome -> hex colour (converted to RGBA)
    edge_colors      : dict — statement type -> hex colour (converted to RGBA)
    node_size        : int  — vertex size (default 28)
    font_size        : int  — vertex label font size (default 8)
    border_width     : float — not used by graph-tool (vertex_color controls border)
    edge_width       : float — edge pen width (default 2.0)
    background_color : str  — hex colour for image background (default '#FFFFFF')

    Returns
    -------
    str — output file path.
    """
    import graph_tool.all as gt

    if G is None:
        G = _build_graph(stmts, reg)

    # Resolve style kwargs
    chrom_colors_hex = {**DEFAULT_CHROM_BG_COLORS, **(kwargs.get('chrom_colors') or {})}
    edge_colors_hex = {**DEFAULT_EDGE_COLORS, **(kwargs.get('edge_colors') or {})}
    node_size_val = kwargs.get('node_size', 28)
    font_size_val = kwargs.get('font_size', 8)
    edge_width_val = kwargs.get('edge_width', 2.0)
    background_color = kwargs.get('background_color', '#FFFFFF')

    # Build RGBA dicts from hex colours
    rgba_node = {}
    for chrom, hex_c in chrom_colors_hex.items():
        rgba_node[chrom] = _hex_to_rgba(hex_c)

    rgba_edge = {}
    for stype, hex_c in edge_colors_hex.items():
        rgba_edge[stype] = _hex_to_rgba(hex_c, alpha=0.8)

    # Vertex fill colour from chromosome
    vfill = G.new_vertex_property('vector<double>')
    for v in G.vertices():
        chrom = G.vp['chromosome'][v]
        vfill[v] = rgba_node.get(chrom, rgba_node.get('unknown', DEFAULT_RGBA_NODE['unknown']))

    # Edge colour and dash
    ecolor = G.new_edge_property('vector<double>')
    edash = G.new_edge_property('vector<double>')
    for e in G.edges():
        st = G.ep['stmt_type'][e]
        ecolor[e] = rgba_edge.get(st, _FALLBACK_RGBA_EDGE)
        edash[e] = [0.02, 0.01] if st == 'Complex' else []

    bg_rgba = _hex_to_rgba(background_color)

    pos = gt.sfdp_layout(G)
    gt.graph_draw(
        G, pos=pos,
        vertex_fill_color=vfill,
        vertex_color=[0.2, 0.2, 0.2, 0.6],
        vertex_size=node_size_val,
        vertex_text=G.vp['name'],
        vertex_text_color=[1, 1, 1, 1],
        vertex_font_size=font_size_val,
        vertex_font_weight=1,
        edge_color=ecolor,
        edge_pen_width=edge_width_val,
        edge_dash_style=edash,
        edge_marker_size=12,
        output_size=size,
        output=output,
        bg_color=bg_rgba,
    )
    print(f'Static PNG \u2192 {output}')
    return output
