"""
visualization.py
----------------
Interactive (ipycytoscape) and static (graph-tool) network visualization
for the INDRA interaction store.

Public API
----------
interactive_network(stmts, reg, G=None, highlight=None, **kwargs) -> ipywidgets.VBox
save_static_png(stmts, reg, G=None, output='interaction_network.png', size=(1800,1300), **kwargs) -> str

Mapping kwargs
--------------
fill_by          : str   — registry field for node fill colour (default 'chromosome').
                   Multi-valued fields (e.g. 'groups', 'contexts') render as pie charts.
border_by        : str   — registry field for node border colour (default 'groups').
                   Must be single-valued; multi-valued fields raise ValueError.
fill_colors      : dict  — {value: hex_color} generic overrides for fill_by colours
border_colors    : dict  — {value: hex_color} generic overrides for border_by colours
chromosome_colors: dict  — {value: hex_color} for chromosome field (applied when
                   fill_by or border_by is 'chromosome')
group_colors     : dict  — {group_name: hex_color} for groups field (applied when
                   fill_by or border_by is 'groups')
context_colors   : dict  — {context_tag: hex_color} for contexts field (applied when
                   fill_by or border_by is 'contexts')

Style kwargs (all optional, applied to both functions where applicable)
-----------------------------------------------------------------------
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

import math
import base64

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


# ── field resolution helpers ─────────────────────────────────────────────────

def _resolve_field(name, reg, field):
    """Extract field values from registry for a gene. Returns list of strings.

    Supports dotted paths (e.g. 'analysis_origin.source') and list fields
    (e.g. 'groups'). Always returns a list; single-valued fields return a
    one-element list.
    """
    entry = reg.get(name, {})
    for part in field.split('.'):
        if isinstance(entry, dict):
            entry = entry.get(part)
        else:
            return ['unknown']
    if entry is None:
        return ['unknown']
    if isinstance(entry, list):
        return entry if entry else ['unknown']
    return [str(entry)]


def _resolve_border_field(name, reg, field):
    """Extract a single border value from registry for a gene.

    For multi-valued fields like 'groups', returns the primary module
    (first matching from MODULE_PRIORITY) to produce a single value.
    For single-valued fields, delegates to _resolve_field.
    """
    if field == 'groups':
        return _primary_module(name, reg)
    vals = _resolve_field(name, reg, field)
    return vals[0]


def _is_multi_valued(reg, field):
    """Return True if *any* gene in the registry has a list for *field*."""
    for name in reg:
        vals = _resolve_field(name, reg, field)
        if len(vals) > 1:
            return True
    # Also check the raw type for the first non-None entry
    for entry in reg.values():
        obj = entry
        for part in field.split('.'):
            if isinstance(obj, dict):
                obj = obj.get(part)
            else:
                obj = None
                break
        if isinstance(obj, list):
            return True
    return False


def _auto_colors(values, palette_name='Set2'):
    """Assign colours from a matplotlib categorical palette to distinct values.

    Returns {value: '#RRGGBB'}.
    """
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap(palette_name)
    sorted_vals = sorted(set(values))
    n = max(len(sorted_vals), 1)
    colors = {}
    for i, val in enumerate(sorted_vals):
        rgba = cmap(i / n if n > 1 else 0.5)
        colors[val] = '#{:02x}{:02x}{:02x}'.format(
            int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255))
    return colors


def _pie_svg(colors, size=64):
    """Generate a data:image/svg+xml URI with equal-sized pie slices.

    *colors* is a list of hex colour strings. Single colour → solid circle.
    Empty → neutral gray circle.
    """
    if not colors:
        colors = ['#DADDE1']
    if len(colors) == 1:
        svg = (
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{size}" height="{size}" viewBox="0 0 {size} {size}">'
            f'<circle cx="{size//2}" cy="{size//2}" r="{size//2}" '
            f'fill="{colors[0]}"/></svg>'
        )
    else:
        cx = cy = size / 2
        r = size / 2
        n = len(colors)
        parts = []
        parts.append(
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{size}" height="{size}" viewBox="0 0 {size} {size}">'
        )
        for i, c in enumerate(colors):
            start_angle = 2 * math.pi * i / n - math.pi / 2
            end_angle = 2 * math.pi * (i + 1) / n - math.pi / 2
            x1 = cx + r * math.cos(start_angle)
            y1 = cy + r * math.sin(start_angle)
            x2 = cx + r * math.cos(end_angle)
            y2 = cy + r * math.sin(end_angle)
            large_arc = 1 if (end_angle - start_angle) > math.pi else 0
            parts.append(
                f'<path d="M{cx},{cy} L{x1:.2f},{y1:.2f} '
                f'A{r},{r} 0 {large_arc} 1 {x2:.2f},{y2:.2f} Z" '
                f'fill="{c}"/>'
            )
        parts.append('</svg>')
        svg = ''.join(parts)
    encoded = base64.b64encode(svg.encode('utf-8')).decode('ascii')
    return f'data:image/svg+xml;base64,{encoded}'


_GRAY = '#DADDE1'


def _resolve_colors(reg, field, user_colors, builtin_defaults, *, border=False):
    """Build a complete value→colour mapping for a field.

    When builtin_defaults or user_colors are provided, only explicitly
    listed values get colours — all others default to gray.  When neither
    is provided, colours are auto-generated from a categorical palette.

    When *border* is True, uses _resolve_border_field (single-valued) to
    collect the values that actually appear.
    """
    # Collect all distinct values across registry
    all_vals = set()
    for name in reg:
        if border:
            all_vals.add(_resolve_border_field(name, reg, field))
        else:
            all_vals.update(_resolve_field(name, reg, field))
    all_vals.add('unknown')

    has_explicit = bool(builtin_defaults) or bool(user_colors)
    if has_explicit:
        # Only explicitly listed values get colour; rest gray
        colors = {v: _GRAY for v in all_vals}
    else:
        # No guidance at all — auto-generate a palette
        colors = _auto_colors(all_vals)

    colors.update(builtin_defaults)
    if user_colors:
        colors.update(user_colors)
    return colors


# Map of field names → their built-in default color dicts
_FIELD_BUILTIN_COLORS = {
    'chromosome': DEFAULT_CHROM_BG_COLORS,
    'groups':     DEFAULT_MODULE_BORDER_COLORS,
}


def _builtin_for_field(field, chromosome_colors=None, group_colors=None,
                       context_colors=None):
    """Return the builtin color dict for a given field.

    When a field-specific dict is provided (e.g. chromosome_colors),
    it *replaces* the hard-coded defaults — only values explicitly
    listed get colours; everything else falls through to gray.
    When no field-specific dict is provided, hard-coded defaults apply.
    """
    field_specific = {
        'chromosome': chromosome_colors,
        'groups':     group_colors,
        'contexts':   context_colors,
    }.get(field)

    if field_specific is not None:
        return dict(field_specific)
    return dict(_FIELD_BUILTIN_COLORS.get(field, {}))


# ── private helpers ──────────────────────────────────────────────────────────

def _filter_stmts(stmts, genes=None):
    """Filter statements to those involving any gene in *genes*.

    Parameters
    ----------
    stmts : list[Statement]
        INDRA statement objects.
    genes : list[str] | set[str] | None
        Gene names to keep. If None, return all statements unchanged.

    Returns filtered list of statements.
    """
    if genes is None:
        return stmts
    gene_set = set(genes)
    return [s for s in stmts
            if any(a is not None and a.name in gene_set
                   for a in s.agent_list())]


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
                    ev_text = (s.evidence[0].text or '')[:200] if s.evidence else ''
                    ctx = s.evidence[0].annotations.get('context', '') if s.evidence else ''
                    edge_evidence[key] = (ev_text, ctx)
        elif len(agents) >= 2:
            key = (agents[0].name, agents[-1].name, stype)
            ev_text = (s.evidence[0].text or '')[:200] if s.evidence else ''
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


def _cytoscape_json(G, reg, edge_evidence,
                    fill_by, border_by, fill_colors, border_colors,
                    fill_multi):
    """Convert graph-tool Graph to {nodes, edges} dict for ipycytoscape.

    When *fill_multi* is True, fill is rendered via pie SVG background images
    and each node gets a 'bg_image' data field.
    """
    nodes = []
    edges = []

    for v in G.vertices():
        name = G.vp['name'][v]
        chrom = G.vp['chromosome'][v]
        rescue = reg.get(name, {}).get('rescue_logic', 'none')
        groups = ', '.join(reg.get(name, {}).get('groups', []))
        info = f'{name}\nChromosome: {chrom}\nGroups: {groups}\nRescue: {rescue}'

        fill_vals = _resolve_field(name, reg, fill_by)
        border_val = _resolve_border_field(name, reg, border_by)

        data = {
            'id': name,
            'label': name,
            'border_val': border_val,
            '_info': info,
        }

        if fill_multi:
            # Only include slices for values that have an explicit colour
            slice_colors = [fill_colors[fv] for fv in fill_vals
                            if fill_colors.get(fv, _GRAY) != _GRAY]
            data['bg_image'] = _pie_svg(slice_colors)
        else:
            data['fill_val'] = fill_vals[0]

        nodes.append({'data': data})

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


def _build_style(fill_colors, border_colors, edge_colors,
                 node_size, node_shape, font_size, border_width, edge_width,
                 text_color, font_weight, text_outline, fill_multi):
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
        'background-color': fill_colors.get('unknown', '#DADDE1'),
        'border-width': border_width,
        'border-color': '#999',
    }
    if fill_multi:
        node_style['background-image'] = 'data(bg_image)'
        node_style['background-fit'] = 'cover'
        node_style['background-color'] = '#DADDE1'
    if text_outline:
        node_style['text-outline-color'] = '#333'
        node_style['text-outline-width'] = 1
    else:
        node_style['text-outline-width'] = 0

    style = [
        {'selector': 'node', 'style': node_style},
    ]

    # Fill colouring (only when single-valued — pie charts handled via bg_image)
    if not fill_multi:
        for val, color in fill_colors.items():
            if val == 'unknown':
                continue
            style.append({
                'selector': f'node[fill_val = "{val}"]',
                'style': {'background-color': color},
            })

    # Border colours
    for val, color in border_colors.items():
        style.append({
            'selector': f'node[border_val = "{val}"]',
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


def _build_legend_html(fill_by, border_by, fill_colors, border_colors,
                       edge_colors, fill_multi):
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

    # Node fill
    fill_label = fill_by.replace('.', ' \u203a ')
    pie_note = ' (pie chart)' if fill_multi else ''
    parts.append(f'<b>Node fill \u2014 {fill_label}{pie_note}</b><br>')
    for val in sorted(fill_colors):
        parts.append(swatch.format(color=fill_colors[val], border='#666') + f'{val}<br>')

    # Node border
    border_label = border_by.replace('.', ' \u203a ')
    parts.append(f'<br><b>Node border \u2014 {border_label}</b><br>')
    for val in sorted(border_colors):
        parts.append(swatch.format(color='#ccc', border=border_colors[val]) + f'{val}<br>')

    # Edges
    parts.append('<br><b>Edges \u2014 Interaction type</b><br>')
    seen = set()
    for stype in ['Activation', 'Inhibition', 'Phosphorylation',
                   'IncreaseAmount', 'DecreaseAmount', 'Complex']:
        color = edge_colors.get(stype, '#AAAAAA')
        dash = 'dashed' if stype == 'Complex' else 'solid'
        arrow = arrow_label.get(stype, 'triangle')
        label = f'{stype} \u2014 {arrow}'
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

def interactive_network(stmts, reg, G=None, highlight=None, genes=None,
                        fill_by='chromosome', border_by='groups',
                        fill_colors=None, border_colors=None,
                        chromosome_colors=None, group_colors=None,
                        context_colors=None, **kwargs):
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
    genes : list[str], optional
        Restrict graph to statements involving these genes.
    fill_by : str
        Registry field for node fill colour (default 'chromosome').
        Multi-valued fields (e.g. 'groups', 'contexts') render as pie charts.
    border_by : str
        Registry field for node border colour (default 'groups').
        Must be single-valued; multi-valued fields raise ValueError.
    fill_colors : dict, optional
        {value: hex_color} generic overrides for fill_by colours.
    border_colors : dict, optional
        {value: hex_color} generic overrides for border_by colours.
    chromosome_colors : dict, optional
        {value: hex_color} for chromosome values (e.g. 'auto', 'X', 'Y').
        Applied whenever fill_by or border_by is 'chromosome'.
    group_colors : dict, optional
        {group_name: hex_color} for gene groups.
        Applied whenever fill_by or border_by is 'groups'.
    context_colors : dict, optional
        {context_tag: hex_color} for context tags.
        Applied whenever fill_by or border_by is 'contexts'.

    Style kwargs (all optional)
    ---------------------------
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

    # Validate border_by: multi-valued fields (other than 'groups', which
    # uses primary-module logic) cannot be used for borders.
    if border_by != 'groups' and _is_multi_valued(reg, border_by):
        raise ValueError(
            f"border_by='{border_by}' is a multi-valued field and cannot be "
            f"used for borders. Use it for fill_by instead, or choose a "
            f"single-valued field like 'chromosome' or 'rescue_logic'."
        )

    fill_multi = _is_multi_valued(reg, fill_by)

    stmts = _filter_stmts(stmts, genes)
    if G is None or genes is not None:
        G = _build_graph(stmts, reg)

    # Backward compat: chrom_colors → chromosome_colors, module_colors → group_colors
    if kwargs.get('chrom_colors') and chromosome_colors is None:
        chromosome_colors = kwargs['chrom_colors']
    if kwargs.get('module_colors') and group_colors is None:
        group_colors = kwargs['module_colors']

    # Resolve colour mappings — field-specific dicts are layered into builtins
    _fc = dict(chromosome_colors=chromosome_colors, group_colors=group_colors,
               context_colors=context_colors)
    fill_builtin = _builtin_for_field(fill_by, **_fc)
    border_builtin = _builtin_for_field(border_by, **_fc)
    resolved_fill = _resolve_colors(reg, fill_by, fill_colors, fill_builtin)
    resolved_border = _resolve_colors(reg, border_by, border_colors, border_builtin, border=True)

    edge_colors_resolved = {**DEFAULT_EDGE_COLORS, **(kwargs.get('edge_colors') or {})}
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
    cyto_json = _cytoscape_json(G, reg, edge_evidence,
                                fill_by, border_by,
                                resolved_fill, resolved_border,
                                fill_multi)

    _box_style = (
        'padding:8px; font-size:13px; color:#333; '
        'background:#f8f8f8; border:1px solid #ddd; border-radius:4px; '
        'min-height:40px; font-family:monospace; white-space:pre-wrap;'
    )
    legend_html = _build_legend_html(fill_by, border_by,
                                     resolved_fill, resolved_border,
                                     edge_colors_resolved, fill_multi)
    info_box = widgets.HTML(
        value=f'<div style="{_box_style}">Click a node or edge for details, or press Show Legend.</div>',
    )

    cyto = ipycytoscape.CytoscapeWidget()
    cyto.graph.add_graph_from_json(cyto_json, directed=True)
    cyto.set_style(_build_style(
        resolved_fill, resolved_border, edge_colors_resolved,
        node_size, node_shape, font_size, border_width, edge_width,
        text_color, font_weight, text_outline, fill_multi,
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
                    size=(1800, 1300), genes=None,
                    fill_by='chromosome', border_by='groups',
                    fill_colors=None, border_colors=None,
                    chromosome_colors=None, group_colors=None,
                    context_colors=None, **kwargs):
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
    genes : list[str], optional
        Restrict graph to statements involving these genes.
    fill_by : str
        Registry field for node fill colour (default 'chromosome').
        Multi-valued fields use the first value with a warning.
    border_by : str
        Registry field for node border (outline) colour (default 'groups').
        Must be single-valued; multi-valued fields raise ValueError.
    fill_colors : dict, optional
        {value: hex_color} generic overrides for fill_by colours.
    border_colors : dict, optional
        {value: hex_color} generic overrides for border_by colours.
    chromosome_colors : dict, optional
        {value: hex_color} for chromosome values. Applied when fill_by or
        border_by is 'chromosome'.
    group_colors : dict, optional
        {group_name: hex_color} for gene groups. Applied when fill_by or
        border_by is 'groups'.
    context_colors : dict, optional
        {context_tag: hex_color} for context tags. Applied when fill_by or
        border_by is 'contexts'.

    Style kwargs (all optional)
    ---------------------------
    edge_colors      : dict — statement type -> hex colour (converted to RGBA)
    node_size        : int  — vertex size (default 28)
    font_size        : int  — vertex label font size (default 8)
    edge_width       : float — edge pen width (default 2.0)
    background_color : str  — hex colour for image background (default '#FFFFFF')

    Returns
    -------
    str — output file path.
    """
    import graph_tool.all as gt
    import warnings

    # Validate border_by: multi-valued fields (other than 'groups', which
    # uses primary-module logic) cannot be used for borders.
    if border_by != 'groups' and _is_multi_valued(reg, border_by):
        raise ValueError(
            f"border_by='{border_by}' is a multi-valued field and cannot be "
            f"used for borders. Use it for fill_by instead, or choose a "
            f"single-valued field like 'chromosome' or 'rescue_logic'."
        )

    fill_multi = _is_multi_valued(reg, fill_by)
    if fill_multi:
        warnings.warn(
            f"fill_by='{fill_by}' is multi-valued; graph-tool does not support "
            f"pie charts — using the first value for each gene.",
            stacklevel=2,
        )

    stmts = _filter_stmts(stmts, genes)
    if G is None or genes is not None:
        G = _build_graph(stmts, reg)

    # Backward compat: chrom_colors → chromosome_colors, module_colors → group_colors
    if kwargs.get('chrom_colors') and chromosome_colors is None:
        chromosome_colors = kwargs['chrom_colors']
    if kwargs.get('module_colors') and group_colors is None:
        group_colors = kwargs['module_colors']

    # Resolve colour mappings — field-specific dicts are layered into builtins
    _fc = dict(chromosome_colors=chromosome_colors, group_colors=group_colors,
               context_colors=context_colors)
    fill_builtin = _builtin_for_field(fill_by, **_fc)
    border_builtin = _builtin_for_field(border_by, **_fc)
    resolved_fill = _resolve_colors(reg, fill_by, fill_colors, fill_builtin)
    resolved_border = _resolve_colors(reg, border_by, border_colors, border_builtin, border=True)

    edge_colors_hex = {**DEFAULT_EDGE_COLORS, **(kwargs.get('edge_colors') or {})}
    node_size_val = kwargs.get('node_size', 28)
    font_size_val = kwargs.get('font_size', 8)
    edge_width_val = kwargs.get('edge_width', 2.0)
    background_color = kwargs.get('background_color', '#FFFFFF')

    # Build RGBA dicts from hex colours
    rgba_fill = {val: _hex_to_rgba(hx) for val, hx in resolved_fill.items()}
    rgba_border = {val: _hex_to_rgba(hx) for val, hx in resolved_border.items()}
    rgba_edge = {st: _hex_to_rgba(hx, alpha=0.8) for st, hx in edge_colors_hex.items()}

    # Vertex fill colour
    vfill = G.new_vertex_property('vector<double>')
    for v in G.vertices():
        name = G.vp['name'][v]
        fill_vals = _resolve_field(name, reg, fill_by)
        val = fill_vals[0]
        vfill[v] = rgba_fill.get(val, rgba_fill.get('unknown', DEFAULT_RGBA_NODE['unknown']))

    # Vertex border (outline) colour
    vcolor = G.new_vertex_property('vector<double>')
    for v in G.vertices():
        name = G.vp['name'][v]
        val = _resolve_border_field(name, reg, border_by)
        vcolor[v] = rgba_border.get(val, [0.2, 0.2, 0.2, 0.6])

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
        vertex_color=vcolor,
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
