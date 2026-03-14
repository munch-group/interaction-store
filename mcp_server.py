"""
mcp_server.py
-------------
Local MCP server for the INDRA interaction store.

Exposes the statement store and gene registry as Claude Code tools,
so the agent can read/write interactions and metadata directly
without file I/O in the conversation loop.

Usage:
    python mcp_server.py            # stdio transport (Claude Code default)
    python mcp_server.py --port 8765  # SSE transport for debugging

Claude Code config (.mcp.json in project root):
    {
      "mcpServers": {
        "interaction-store": {
          "command": "pixi",
          "args": ["run", "python", "mcp_server.py"],
          "cwd": "${workspaceFolder}"
        }
      }
    }

If not using pixi, replace with:
    "command": "python",
    "args": ["mcp_server.py"]
"""

import argparse
import json
import pathlib
import sys
from datetime import date
from typing import Optional

from fastmcp import FastMCP

# ── paths ─────────────────────────────────────────────────────────────────────
BASE          = pathlib.Path(__file__).parent
STORE_PATH    = BASE / 'statements.json'
REGISTRY_PATH = BASE / 'gene_registry.json'

# ── lazy INDRA imports (only needed at call time) ─────────────────────────────
def _indra():
    from indra.statements import (
        Agent, Evidence,
        Activation, Inhibition,
        Phosphorylation, Dephosphorylation,
        Ubiquitination, Deubiquitination,
        Acetylation, Deacetylation,
        Methylation, Demethylation,
        IncreaseAmount, DecreaseAmount,
        Complex, Translocation,
        stmts_to_json, stmts_from_json,
    )
    from indra.assemblers.indranet import IndraNetAssembler
    return locals()


# ── store I/O ─────────────────────────────────────────────────────────────────
def _load_store():
    from indra.statements import stmts_from_json
    if not STORE_PATH.exists():
        return []
    with open(STORE_PATH) as f:
        return stmts_from_json(json.load(f))


def _save_store(stmts):
    from indra.statements import stmts_to_json
    with open(STORE_PATH, 'w') as f:
        json.dump(stmts_to_json(stmts), f, indent=2)


def _load_registry() -> dict:
    if not REGISTRY_PATH.exists():
        return {}
    with open(REGISTRY_PATH) as f:
        return json.load(f)


def _save_registry(reg: dict):
    with open(REGISTRY_PATH, 'w') as f:
        json.dump(reg, f, indent=2, sort_keys=True)


def _make_agent(name: str, hgnc: str = None, up: str = None):
    from indra.statements import Agent
    db_refs = {}
    if hgnc: db_refs['HGNC'] = str(hgnc)
    if up:   db_refs['UP']   = up
    return Agent(name, db_refs=db_refs or None)


def _make_evidence(text: str, pmid: str = None, doi: str = None,
                   context: str = None, hypothesis: bool = False,
                   direct: bool = True, species: str = 'Homo_sapiens'):
    from indra.statements import Evidence
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
            'species':    species,
        },
        epistemics={'hypothesis': hypothesis},
    )


# ── genome coordinates ─────────────────────────────────────────────────────────

def _lookup_coords(gene_name: str) -> dict:
    """Look up genome coordinates. Currently supports hg38 via geneinfo;
    additional assemblies can be added here as geneinfo gains support."""
    coords = {}
    try:
        from geneinfo.coords import gene_coords
        hits = gene_coords(gene_name, 'hg38')
        if hits:
            chrom, start, end, _ = hits[0]
            coords['hg38'] = {'chrom': chrom, 'start': start, 'end': end}
    except Exception:
        pass
    return coords


# ── server ────────────────────────────────────────────────────────────────────
mcp = FastMCP(
    name='interaction-store',
    instructions=(
        'Tools for managing the INDRA mechanistic interaction store '
        'and gene attribute registry for exploratory research. '
        'Refer to CLAUDE.md in the project root for statement types, '
        'evidence conventions, and gene group definitions.'
    ),
)


# ════════════════════════════════════════════════════════════════════════════
# STATEMENT TOOLS
# ════════════════════════════════════════════════════════════════════════════

@mcp.tool
def add_statement(
    stmt_type: str,
    subject: str,
    object: str,
    evidence_text: str,
    pmid: Optional[str] = None,
    doi: Optional[str] = None,
    context: Optional[str] = None,
    hypothesis: bool = False,
    direct: bool = True,
    residue: Optional[str] = None,
    position: Optional[str] = None,
    subject_hgnc: Optional[str] = None,
    object_hgnc: Optional[str] = None,
    species: str = 'Homo_sapiens',
) -> str:
    """
    Add a new INDRA statement to the persistent store.

    stmt_type must be one of the types defined in CLAUDE.md §1:
      Activity:     Activation, Inhibition
      Amount:       IncreaseAmount, DecreaseAmount
      Modification: Phosphorylation, Dephosphorylation, Ubiquitination,
                    Deubiquitination, Acetylation, Deacetylation,
                    Methylation, Demethylation, (+ all others in §1.2)
      Other:        Complex, Translocation, Gef, Gap

    For Complex, subject and object are the two primary members.
    For modifications, subject=enzyme, object=substrate.
    species defaults to 'Homo_sapiens'. Use e.g. 'Mus_musculus',
    'Pan_troglodytes' for non-human evidence.
    """
    import indra.statements as IS

    cls = getattr(IS, stmt_type, None)
    if cls is None:
        return f'ERROR: unknown statement type "{stmt_type}". See CLAUDE.md §1.'

    subj = _make_agent(subject, hgnc=subject_hgnc)
    obj  = _make_agent(object,  hgnc=object_hgnc)
    ev   = _make_evidence(evidence_text, pmid=pmid, doi=doi,
                          context=context, hypothesis=hypothesis,
                          direct=direct, species=species)

    # Build statement — handle type families
    mod_types = {
        'Phosphorylation', 'Dephosphorylation', 'Ubiquitination',
        'Deubiquitination', 'Acetylation', 'Deacetylation',
        'Methylation', 'Demethylation', 'Sumoylation', 'Desumoylation',
        'Hydroxylation', 'Dehydroxylation', 'Ribosylation', 'Deribosylation',
        'Glycosylation', 'Deglycosylation', 'Farnesylation', 'Defarnesylation',
        'Geranylgeranylation', 'Degeranylgeranylation',
        'Palmitoylation', 'Depalmitoylation',
        'Myristoylation', 'Demyristoylation',
    }
    if stmt_type in mod_types:
        stmt = cls(subj, obj, residue=residue, position=position, evidence=[ev])
    elif stmt_type == 'Complex':
        stmt = IS.Complex([subj, obj], evidence=[ev])
    else:
        stmt = cls(subj, obj, evidence=[ev])

    stmts = _load_store()
    stmts.append(stmt)
    _save_store(stmts)

    ref_str = ''
    if pmid: ref_str += f' PMID:{pmid}'
    if doi:  ref_str += f' DOI:{doi}'
    hyp_str = ' [HYPOTHESIS]' if hypothesis else ''
    dir_str = ' [indirect]' if not direct else ''
    return (
        f'Added {stmt_type}({subject}, {object}){hyp_str}{dir_str}{ref_str}\n'
        f'Store now has {len(stmts)} statements.'
    )


@mcp.tool
def query_statements(
    gene: Optional[str] = None,
    stmt_type: Optional[str] = None,
    context: Optional[str] = None,
    hypothesis_only: bool = False,
) -> str:
    """
    Query the statement store.

    Filter by gene name (subject or object), statement type,
    context annotation, or show only hypothetical statements.
    Returns a human-readable summary with evidence text and references.
    """
    stmts = _load_store()

    if gene:
        stmts = [s for s in stmts
                 if any(a and a.name == gene for a in s.agent_list())]
    if stmt_type:
        stmts = [s for s in stmts if type(s).__name__ == stmt_type]
    if context:
        stmts = [s for s in stmts
                 if any(context.lower() in
                        (e.annotations.get('context') or '').lower()
                        for e in s.evidence)]
    if hypothesis_only:
        stmts = [s for s in stmts
                 if any((e.epistemics or {}).get('hypothesis')
                        for e in s.evidence)]

    if not stmts:
        return 'No statements match the query.'

    lines = [f'{len(stmts)} statement(s) found:\n']
    for s in stmts:
        agents = [a.name for a in s.agent_list() if a]
        ev = s.evidence[0]
        refs = []
        for e in s.evidence:
            tr = e.text_refs or {}
            if tr.get('DOI'):  refs.append(f"DOI:{tr['DOI']}")
            elif tr.get('PMID'): refs.append(f"PMID:{tr['PMID']}")
        hyp = any((e.epistemics or {}).get('hypothesis') for e in s.evidence)
        ctx = ev.annotations.get('context', '')
        drct = ev.annotations.get('directness', '')
        lines.append(
            f'[{type(s).__name__}] {" → ".join(agents)}\n'
            f'  context:    {ctx}\n'
            f'  directness: {drct}{"  [HYPOTHESIS]" if hyp else ""}\n'
            f'  evidence:   {ev.text[:120]}\n'
            f'  refs:       {", ".join(refs) if refs else "none"}\n'
        )
    return '\n'.join(lines)


@mcp.tool
def promote_to_literature(
    subject: str,
    object: str,
    stmt_type: str,
    evidence_text: str,
    pmid: Optional[str] = None,
    doi: Optional[str] = None,
    context: Optional[str] = None,
) -> str:
    """
    Append a literature-backed Evidence to an existing statement,
    promoting it from hypothesis to supported.

    Matches the first statement of the given type with matching
    subject and object names. Keeps all existing evidence intact.
    """
    stmts = _load_store()
    for s in stmts:
        agents = [a.name for a in s.agent_list() if a]
        if (type(s).__name__ == stmt_type
                and len(agents) >= 2
                and agents[0] == subject
                and agents[-1] == object):
            ev = _make_evidence(
                evidence_text, pmid=pmid, doi=doi,
                context=context, hypothesis=False,
            )
            s.evidence.append(ev)
            _save_store(stmts)
            ref = doi or pmid or 'no ref'
            return (
                f'Appended literature evidence to '
                f'{stmt_type}({subject}, {object}).\n'
                f'Reference: {ref}\n'
                f'Statement now has {len(s.evidence)} evidence object(s).'
            )
    return (
        f'No matching {stmt_type}({subject}, {object}) found. '
        f'Use add_statement to create it first.'
    )


@mcp.tool
def list_statements_summary() -> str:
    """
    Return a compact summary of the full statement store:
    counts by type, hypothesis count, and genes covered.
    """
    stmts = _load_store()
    if not stmts:
        return 'Store is empty.'

    by_type: dict[str, int] = {}
    genes: set[str] = set()
    hyp_count = 0
    ref_count = 0

    for s in stmts:
        t = type(s).__name__
        by_type[t] = by_type.get(t, 0) + 1
        for a in s.agent_list():
            if a: genes.add(a.name)
        if any((e.epistemics or {}).get('hypothesis') for e in s.evidence):
            hyp_count += 1
        if any((e.text_refs or {}) for e in s.evidence):
            ref_count += 1

    lines = [
        f'Total statements: {len(stmts)}',
        f'  with DOI/PMID reference: {ref_count}',
        f'  flagged as hypothesis:   {hyp_count}',
        '',
        'By type:',
    ]
    for t, n in sorted(by_type.items()):
        lines.append(f'  {t:25s} {n}')
    lines += ['', f'Genes covered ({len(genes)}):',
              '  ' + ', '.join(sorted(genes))]
    return '\n'.join(lines)


@mcp.tool
def enrich_from_indra_db(
    genes: list[str],
    max_age_days: Optional[int] = None,
) -> str:
    """
    Query the public INDRA DB for known statements involving the given
    genes. Uses a persistent per-gene cache (indra_db_cache.json) so
    repeated queries are instant.

    Parameters
    ----------
    genes : list[str]
        Gene names to query (one API call per gene).
    max_age_days : int, optional
        Re-fetch cached entries older than this. Default: never expire.

    Returns a summary — does NOT auto-persist to the statement store.
    """
    from indra_cache import get_statements_batch, cache_summary

    try:
        stmts = get_statements_batch(
            genes, ev_limit=5, max_age_days=max_age_days,
            sleep_between=1.0, verbose=False,
        )
    except Exception as e:
        return f'INDRA DB query failed: {e}'

    if not stmts:
        return 'No statements found in INDRA DB for the given genes.'

    by_type: dict[str, int] = {}
    for s in stmts:
        t = type(s).__name__
        by_type[t] = by_type.get(t, 0) + 1

    lines = [
        f'Found {len(stmts)} statements in INDRA DB for {genes}',
        '',
        'Breakdown by type:',
    ]
    for t, n in sorted(by_type.items()):
        lines.append(f'  {t:25s} {n}')
    lines += [
        '',
        'Top 20 statements:',
    ]
    for s in stmts[:20]:
        agents = [a.name for a in s.agent_list() if a]
        pmids = list({e.pmid for e in s.evidence if e.pmid})[:3]
        lines.append(
            f'  [{type(s).__name__:20s}] {" → ".join(agents):40s}'
            f'  PMIDs: {pmids}'
        )
    if len(stmts) > 20:
        lines.append(f'  ... and {len(stmts) - 20} more')
    lines.append(
        '\nTo persist selected statements, use add_statement for '
        'each interaction you want to keep.'
    )
    lines.append(f'\n{cache_summary()}')
    return '\n'.join(lines)


# ════════════════════════════════════════════════════════════════════════════
# REGISTRY TOOLS
# ════════════════════════════════════════════════════════════════════════════

@mcp.tool
def add_gene_to_registry(
    name: str,
    chromosome: str,
    groups: Optional[list[str]] = None,
    analysis_source: Optional[str] = None,
    analysis_name: Optional[str] = None,
    analysis_note: Optional[str] = None,
    rescue_logic: Optional[str] = None,
    contexts: Optional[list[str]] = None,
    haplogroup_effect: Optional[str] = None,
    notes: Optional[str] = None,
    pmid: Optional[str] = None,
    doi: Optional[str] = None,
    reference_note: Optional[str] = None,
) -> str:
    """
    Add or update a gene entry in the registry.

    chromosome: 'auto', 'X', 'Y', or 'mito'
    rescue_logic: 'rheostat', 'paralog_backup', 'partial', or 'none'
    analysis_source: e.g. 'IBDmix_NHR', 'GWAS', 'literature', 'manual'
    """
    reg = _load_registry()
    entry = reg.get(name, {})

    entry['chromosome'] = chromosome
    if groups:
        existing = entry.get('groups', [])
        for g in groups:
            if g not in existing:
                existing.append(g)
        entry['groups'] = existing
    if analysis_source or analysis_name or analysis_note:
        entry['analysis_origin'] = {
            k: v for k, v in {
                'source':   analysis_source,
                'analysis': analysis_name,
                'note':     analysis_note,
            }.items() if v
        }
    if rescue_logic:
        entry['rescue_logic'] = rescue_logic
    if contexts:
        entry['contexts'] = contexts
    if haplogroup_effect:
        entry['haplogroup_effect'] = haplogroup_effect
    if notes:
        entry['notes'] = notes
    # Auto-populate genome coordinates if not already present
    if 'coordinates' not in entry:
        coords = _lookup_coords(name)
        if coords:
            entry['coordinates'] = coords
    if pmid or doi:
        refs = entry.get('references', [])
        refs.append({k: v for k, v in {
            'pmid': pmid, 'doi': doi, 'note': reference_note
        }.items() if v})
        entry['references'] = refs

    reg[name] = entry
    _save_registry(reg)
    return f'Registry updated for {name}.\n{json.dumps(entry, indent=2)}'


@mcp.tool
def query_registry(
    gene: Optional[str] = None,
    group: Optional[str] = None,
    chromosome: Optional[str] = None,
    rescue_only: bool = False,
    analysis: Optional[str] = None,
) -> str:
    """
    Query the gene registry.

    Filter by gene name, group membership, chromosome, rescue
    candidacy, or analysis of origin.
    """
    reg = _load_registry()
    items = list(reg.items())

    if gene:
        items = [(n, a) for n, a in items if n == gene]
    if group:
        items = [(n, a) for n, a in items
                 if group in a.get('groups', [])]
    if chromosome:
        items = [(n, a) for n, a in items
                 if a.get('chromosome') == chromosome]
    if rescue_only:
        items = [(n, a) for n, a in items
                 if a.get('rescue_logic') not in (None, 'none')]
    if analysis:
        items = [(n, a) for n, a in items
                 if analysis in (a.get('analysis_origin') or {}).get('analysis', '')]

    if not items:
        return 'No registry entries match the query.'

    lines = [f'{len(items)} gene(s) found:\n']
    for name, attrs in items:
        groups = ', '.join(attrs.get('groups', [])) or '—'
        chrom  = attrs.get('chromosome', '?')
        rescue = attrs.get('rescue_logic', 'none')
        origin = (attrs.get('analysis_origin') or {}).get('analysis', '—')
        refs   = attrs.get('references', [])
        ref_str = '; '.join(
            r.get('doi') or r.get('pmid', '') for r in refs
        ) or 'none'
        lines.append(
            f'{name} [{chrom}]\n'
            f'  groups:        {groups}\n'
            f'  rescue_logic:  {rescue}\n'
            f'  origin:        {origin}\n'
            f'  refs:          {ref_str}\n'
            f'  notes:         {attrs.get("notes", "")[:100]}\n'
        )
    return '\n'.join(lines)


@mcp.tool
def get_gene_group(group_name: str) -> str:
    """Return all genes belonging to a named group."""
    reg = _load_registry()
    members = [n for n, a in reg.items()
               if group_name in a.get('groups', [])]
    if not members:
        return f'No genes found in group "{group_name}".'
    return f'Group "{group_name}" ({len(members)} genes):\n  ' + \
           ', '.join(sorted(members))


@mcp.tool
def list_registry_summary() -> str:
    """
    Return a full summary of the registry:
    genes by chromosome, all groups and members,
    rescue candidates, and analysis provenance.
    """
    reg = _load_registry()
    if not reg:
        return 'Registry is empty.'

    by_chrom: dict[str, list] = {}
    for name, attrs in reg.items():
        c = attrs.get('chromosome', 'unknown')
        by_chrom.setdefault(c, []).append(name)

    all_groups: dict[str, list] = {}
    for name, attrs in reg.items():
        for g in attrs.get('groups', []):
            all_groups.setdefault(g, []).append(name)

    rescue = [n for n, a in reg.items()
              if a.get('rescue_logic') not in (None, 'none')]

    origins: dict[str, list] = {}
    for name, attrs in reg.items():
        src = (attrs.get('analysis_origin') or {}).get('analysis', '')
        if src:
            origins.setdefault(src, []).append(name)

    lines = [f'Registry: {len(reg)} genes\n']
    for chrom in ('auto', 'X', 'Y', 'mito', 'unknown'):
        genes = sorted(by_chrom.get(chrom, []))
        if genes:
            lines.append(f'{chrom:8s}: {", ".join(genes)}')

    lines.append('\nGroups:')
    for g, members in sorted(all_groups.items()):
        lines.append(f'  [{g}]\n    {", ".join(sorted(members))}')

    lines.append(f'\nRescue candidates: {", ".join(sorted(rescue))}')

    if origins:
        lines.append('\nAnalysis provenance:')
        for analysis, genes in sorted(origins.items()):
            lines.append(f'  {analysis}: {", ".join(sorted(genes))}')

    return '\n'.join(lines)


# ════════════════════════════════════════════════════════════════════════════
# BROWSING TOOLS
# ════════════════════════════════════════════════════════════════════════════

@mcp.tool
def list_contexts() -> str:
    """
    List all context tags used in the statement store with counts.
    """
    from gene_registry import get_all_contexts
    contexts = get_all_contexts(str(STORE_PATH))
    if not contexts:
        return 'No context tags found.'
    lines = [f'{len(contexts)} context tags:\n']
    for ctx, n in sorted(contexts.items(), key=lambda x: -x[1]):
        lines.append(f'  {n:3d}  {ctx}')
    return '\n'.join(lines)


@mcp.tool
def genes_by_context(context: str) -> str:
    """
    List all genes appearing in statements with a given context tag.
    Case-insensitive substring match.
    """
    from gene_registry import genes_by_context as _genes_by_context
    genes = list(_genes_by_context(context, str(STORE_PATH)))
    if not genes:
        return f'No genes found for context "{context}".'
    return (
        f'Context "{context}" — {len(genes)} gene(s):\n'
        f'  {", ".join(genes)}'
    )


@mcp.tool
def gene_interactions(gene: str) -> str:
    """
    List all genes that interact with a specific gene in the store,
    with statement type, context, and references.
    """
    from gene_registry import get_interactors
    hits = get_interactors(gene, str(STORE_PATH))
    if not hits:
        return f'No interactions found for {gene}.'
    lines = [f'{gene} — {len(hits)} interaction(s):\n']
    for h in hits:
        ref_str = ', '.join(h['refs']) if h['refs'] else 'none'
        lines.append(
            f'  [{h["type"]}] {gene} — {h["gene"]}\n'
            f'    context: {h["context"]}\n'
            f'    refs:    {ref_str}'
        )
    return '\n'.join(lines)


@mcp.tool
def gene_info(gene: str) -> str:
    """
    Show all recorded information about a gene: registry entry
    plus all statements involving it.
    """
    from gene_registry import get_gene_info, get_interactors
    info = get_gene_info(gene)
    if info is None:
        reg_section = f'{gene} is not in the registry.'
    else:
        lines = [f'{gene} [{info.get("chromosome", "?")}]']
        if info.get('groups'):
            lines.append(f'  groups:       {", ".join(info["groups"])}')
        if info.get('rescue_logic') and info['rescue_logic'] != 'none':
            lines.append(f'  rescue:       {info["rescue_logic"]}')
        if info.get('contexts'):
            lines.append(f'  expression:   {", ".join(info["contexts"])}')
        if info.get('haplogroup_effect'):
            lines.append(f'  haplogroup:   {info["haplogroup_effect"]}')
        origin = info.get('analysis_origin')
        if origin:
            lines.append(f'  analysis:     {origin.get("analysis", "")} ({origin.get("source", "")})')
            if origin.get('note'):
                lines.append(f'                {origin["note"]}')
        refs = info.get('references', [])
        if refs:
            lines.append(f'  references:')
            for r in refs:
                ref_id = r.get('doi') or r.get('pmid', '?')
                note = r.get('note', '')
                lines.append(f'    {ref_id}  {note}')
        if info.get('notes'):
            lines.append(f'  notes:        {info["notes"]}')
        reg_section = '\n'.join(lines)

    # Interactions
    hits = get_interactors(gene, str(STORE_PATH))
    if hits:
        int_lines = [f'\n{len(hits)} interaction(s):']
        for h in hits:
            ref_str = ', '.join(h['refs']) if h['refs'] else 'none'
            int_lines.append(
                f'  [{h["type"]}] {gene} — {h["gene"]}  ({h["context"]})  refs: {ref_str}'
            )
        interaction_section = '\n'.join(int_lines)
    else:
        interaction_section = '\nNo interactions in the store.'

    return reg_section + interaction_section


@mcp.tool
def list_groups() -> str:
    """List all gene groups in the registry with their members."""
    from gene_registry import get_all_groups
    groups = get_all_groups()
    if not groups:
        return 'No groups found.'
    lines = [f'{len(groups)} group(s):\n']
    for name, members in groups.items():
        lines.append(f'  [{name}] ({len(members)})')
        lines.append(f'    {", ".join(members)}')
    return '\n'.join(lines)


# ════════════════════════════════════════════════════════════════════════════
# GRAPH TOOL
# ════════════════════════════════════════════════════════════════════════════

def _build_gt_graph(stmts, reg):
    """
    Build a graph-tool directed Graph from INDRA statements + registry.

    Returns (G, name_prop, stmt_type_prop) where the property maps are
    also attached as G.vp['name'] and G.ep['stmt_type'].
    """
    from graph_tool import Graph as GTGraph

    G = GTGraph(directed=True)
    name_prop = G.new_vertex_property('string')
    G.vp['name'] = name_prop

    chrom_prop = G.new_vertex_property('string')
    G.vp['chromosome'] = chrom_prop

    stmt_type_prop = G.new_edge_property('string')
    G.ep['stmt_type'] = stmt_type_prop

    # Build node index
    node_idx: dict[str, int] = {}

    def _get_or_add(name: str):
        if name not in node_idx:
            v = G.add_vertex()
            idx = int(v)
            node_idx[name] = idx
            name_prop[v] = name
            info = reg.get(name, {})
            chrom_prop[v] = info.get('chromosome', 'unknown')
        return G.vertex(node_idx[name])

    for s in stmts:
        agents = [a for a in s.agent_list() if a is not None]
        stype = type(s).__name__

        if stype == 'Complex':
            # Undirected: add edges both ways for visibility
            for i, a in enumerate(agents):
                for b in agents[i+1:]:
                    va, vb = _get_or_add(a.name), _get_or_add(b.name)
                    e = G.add_edge(va, vb)
                    stmt_type_prop[e] = 'Complex'
        elif len(agents) >= 2:
            va = _get_or_add(agents[0].name)
            vb = _get_or_add(agents[-1].name)
            e = G.add_edge(va, vb)
            stmt_type_prop[e] = stype

    return G, name_prop, stmt_type_prop


@mcp.tool
def render_network(
    highlight_gene: Optional[str] = None,
    output_path: Optional[str] = None,
    genes: Optional[list[str]] = None,
) -> str:
    """
    Rebuild the interaction network graph from the current store,
    optionally highlighting a specific gene and its neighbours.
    Saves PNG to output_path (default: interaction_network.png).

    If genes is provided, only statements involving those genes are included.
    """
    from graph_tool import Graph as GTGraph
    from graph_tool.draw import graph_draw, sfdp_layout

    stmts = _load_store()
    reg   = _load_registry()
    out   = str(pathlib.Path(output_path or (BASE / 'interaction_network.png')))

    if genes:
        gene_set = set(genes)
        stmts = [s for s in stmts
                 if any(a is not None and a.name in gene_set
                        for a in s.agent_list())]

    G, name_prop, stmt_type_prop = _build_gt_graph(stmts, reg)
    chrom_prop = G.vp['chromosome']

    COLORS = {
        'auto': [0.357, 0.608, 0.835, 1.0],
        'X':    [0.851, 0.373, 0.373, 1.0],
        'Y':    [0.910, 0.627, 0.125, 1.0],
        'unknown': [0.855, 0.867, 0.882, 1.0],
    }
    EDGE_COLORS = {
        'Activation':      [0.153, 0.682, 0.376, 0.8],
        'Inhibition':      [0.753, 0.224, 0.169, 0.8],
        'Phosphorylation': [0.161, 0.502, 0.725, 0.8],
        'IncreaseAmount':  [0.153, 0.682, 0.376, 0.8],
        'DecreaseAmount':  [0.753, 0.224, 0.169, 0.8],
        'Complex':         [0.533, 0.533, 0.533, 0.8],
    }
    DEFAULT_EDGE = [0.667, 0.667, 0.667, 0.8]

    # Vertex fill colour
    vfill = G.new_vertex_property('vector<double>')
    valpha = G.new_vertex_property('double')
    for v in G.vertices():
        c = chrom_prop[v]
        vfill[v] = COLORS.get(c, COLORS['unknown'])
        valpha[v] = 0.92

    # Highlight mode
    if highlight_gene:
        name_to_v = {name_prop[v]: v for v in G.vertices()}
        hv = name_to_v.get(highlight_gene)
        if hv is not None:
            nbrs = set()
            for e in hv.out_edges():
                nbrs.add(int(e.target()))
            for e in hv.in_edges():
                nbrs.add(int(e.source()))
            nbrs.add(int(hv))
            for v in G.vertices():
                if int(v) not in nbrs:
                    vfill[v] = [0.933, 0.933, 0.933, 1.0]
                    valpha[v] = 0.25

    # Edge colour
    ecolor = G.new_edge_property('vector<double>')
    for e in G.edges():
        st = stmt_type_prop[e]
        ecolor[e] = EDGE_COLORS.get(st, DEFAULT_EDGE)

    # Edge dash for Complex
    edash = G.new_edge_property('vector<double>')
    for e in G.edges():
        if stmt_type_prop[e] == 'Complex':
            edash[e] = [0.02, 0.01]
        else:
            edash[e] = []

    # Layout
    pos = sfdp_layout(G)

    graph_draw(
        G, pos=pos,
        vertex_fill_color=vfill,
        vertex_color=[0.2, 0.2, 0.2, 0.6],
        vertex_size=28,
        vertex_text=name_prop,
        vertex_text_color=[1, 1, 1, 1],
        vertex_font_size=8,
        vertex_font_weight=1,  # cairo.FONT_WEIGHT_BOLD
        edge_color=ecolor,
        edge_pen_width=2.0,
        edge_dash_style=edash,
        edge_marker_size=12,
        output_size=(1800, 1300),
        output=out,
    )

    n_nodes = G.num_vertices()
    n_edges = G.num_edges()
    return (
        f'Graph saved → {out}\n'
        f'{n_nodes} nodes, {n_edges} edges.'
    )


# ════════════════════════════════════════════════════════════════════════════
# NL EXTRACTION TOOLS
# ════════════════════════════════════════════════════════════════════════════

PENDING_PATH = BASE / 'pending_extraction.json'


def _format_stmt_preview(idx: int, stmt) -> str:
    """Return a concise numbered preview line for one INDRA statement."""
    from indra.statements import Evidence
    stype  = type(stmt).__name__
    agents = [a.name for a in stmt.agent_list() if a is not None]
    ev0    = stmt.evidence[0] if stmt.evidence else Evidence()
    src    = ev0.source_api or '?'
    pmid   = ev0.pmid or (ev0.text_refs or {}).get('PMID', '')
    ref    = f'  PMID:{pmid}' if pmid else ''

    # Extra detail for PTMs
    detail = ''
    if hasattr(stmt, 'residue') and stmt.residue:
        detail = f' [{stmt.residue}{stmt.position or ""}]'

    belief_str = f'  belief={stmt.belief:.2f}' if hasattr(stmt, 'belief') else ''

    return (
        f'  [{idx}] {stype}({", ".join(agents)}){detail}\n'
        f'       src={src}{ref}{belief_str}\n'
        f'       "{(ev0.text or "")[:120]}"'
    )


@mcp.tool
def extract_from_text(
    text: str,
    reader: str = 'trips',
    pmid: Optional[str] = None,
    context: Optional[str] = None,
    service_host: Optional[str] = None,
) -> str:
    """
    Extract INDRA statements from free text using an NLP reader.

    Sends text to a remote reader service and returns a numbered preview
    of extracted statements. Does NOT persist automatically — use
    persist_extracted() to selectively add statements to the store.

    Parameters
    ----------
    text : str
        One or more sentences describing molecular interactions.
        Works best with precise mechanistic language:
          "GSK3β phosphorylates MAPT at Thr231, reducing its
           affinity for microtubules."
    reader : str
        'trips' (default) — TRIPS/DRUM web service at IHMC.
            Good for single sentences and short passages.
            Endpoint: http://trips.ihmc.us/parser/cgi/drum
        'reach' — REACH/ODIN web service at Arizona.
            Better for full sentences from papers.
            Endpoint: http://agathon.sista.arizona.edu:8080/odinweb/api/text
        'reach_local' — REACH running on localhost:8080.
            Use if you have the Docker container running:
              docker run -d -p 8080:8080 clulab/reach
    pmid : str, optional
        PubMed ID to attach to extracted evidence (if text is from a paper).
    context : str, optional
        Module/context annotation to attach to extracted evidence.
    service_host : str, optional
        Override the default service URL (e.g. for local TRIPS instance).

    Returns
    -------
    Numbered preview of extracted statements. Call persist_extracted()
    with chosen indices to add them to the store.

    Notes
    -----
    Both public endpoints are best-effort research services; availability
    is not guaranteed. For offline/bulk use, run REACH locally via Docker.
    NL readers make grounding and polarity errors — always review before
    persisting.
    """
    reader = reader.strip().lower()

    # ── call the reader ───────────────────────────────────────────────────
    try:
        if reader == 'trips':
            from indra.sources import trips
            kwargs = {}
            if service_host:
                kwargs['service_host'] = service_host
            proc = trips.process_text(text, **kwargs)

        elif reader == 'reach':
            from indra.sources import reach
            url = (service_host or
                   'http://agathon.sista.arizona.edu:8080/odinweb/api/text')
            proc = reach.process_text(text, citation=pmid, url=url)

        elif reader == 'reach_local':
            from indra.sources import reach
            proc = reach.process_text(
                text,
                citation=pmid,
                url='http://localhost:8080/api/text',
            )

        else:
            return (
                f'ERROR: unknown reader "{reader}". '
                'Choose trips, reach, or reach_local.'
            )

    except Exception as exc:
        return (
            f'ERROR: reader call failed ({type(exc).__name__}: {exc}).\n\n'
            'Possible causes:\n'
            '  • Remote service is down or unreachable\n'
            '  • Network access is blocked from this environment\n'
            '  • For reach_local: Docker container not running\n\n'
            'To run REACH locally:\n'
            '  docker run -d -p 8080:8080 clulab/reach\n'
            'Then retry with reader="reach_local".'
        )

    stmts = proc.statements if proc else []

    if not stmts:
        return (
            'Reader returned no statements.\n\n'
            'Tips:\n'
            '  • Use precise mechanistic verbs: phosphorylates, activates,\n'
            '    inhibits, binds, ubiquitinates, recruits\n'
            '  • Name both agents explicitly (avoid pronouns)\n'
            '  • Try one sentence at a time\n'
            '  • TRIPS handles gene/protein names better than pathway names'
        )

    # ── annotate evidence with context/pmid ──────────────────────────────
    from indra.statements import Evidence
    for stmt in stmts:
        for ev in stmt.evidence:
            ev.annotations.setdefault('context', context or 'extracted')
            ev.annotations['extraction_reader'] = reader
            ev.annotations['date'] = str(date.today())
            if pmid and not ev.pmid:
                ev.pmid = pmid
                ev.text_refs = ev.text_refs or {}
                ev.text_refs['PMID'] = str(pmid)

    # ── save pending ──────────────────────────────────────────────────────
    from indra.statements import stmts_to_json
    with open(PENDING_PATH, 'w') as f:
        json.dump(stmts_to_json(stmts), f, indent=2)

    # ── format preview ────────────────────────────────────────────────────
    lines = [
        f'{len(stmts)} statement(s) extracted by {reader} '
        f'(saved to pending_extraction.json):\n'
    ]
    for i, s in enumerate(stmts):
        lines.append(_format_stmt_preview(i, s))
        lines.append('')

    lines += [
        '─' * 60,
        'To persist selected statements:',
        '  persist_extracted(indices=[0, 2, 4])',
        'To persist all:',
        '  persist_extracted()',
        '',
        'Review carefully — NL readers can misidentify agents,',
        'reverse polarity, or miss context-dependent meaning.',
    ]

    return '\n'.join(lines)


@mcp.tool
def persist_extracted(
    indices: Optional[list[int]] = None,
    context: Optional[str] = None,
    hypothesis: bool = False,
) -> str:
    """
    Persist statements from the most recent extract_from_text() call.

    Reads pending_extraction.json written by extract_from_text, adds
    selected statements to the main store, then clears the pending file.

    Parameters
    ----------
    indices : list[int], optional
        Zero-based indices of statements to persist (from the preview).
        If None or empty, all pending statements are persisted.
    context : str, optional
        Override the context annotation on persisted statements.
    hypothesis : bool
        If True, mark all persisted statements as hypothetical
        (useful when promoting NL-extracted statements that need
        manual verification before being treated as confirmed).
    """
    if not PENDING_PATH.exists():
        return (
            'No pending extraction found. '
            'Run extract_from_text() first.'
        )

    from indra.statements import stmts_from_json
    with open(PENDING_PATH) as f:
        pending = stmts_from_json(json.load(f))

    if not pending:
        return 'Pending extraction is empty.'

    # Select
    if indices:
        bad = [i for i in indices if i >= len(pending) or i < 0]
        if bad:
            return (
                f'ERROR: indices {bad} out of range '
                f'(pending has {len(pending)} statements, '
                f'indices 0–{len(pending)-1}).'
            )
        chosen = [pending[i] for i in indices]
    else:
        chosen = pending

    # Apply overrides
    if context or hypothesis:
        for stmt in chosen:
            for ev in stmt.evidence:
                if context:
                    ev.annotations['context'] = context
                if hypothesis:
                    ev.epistemics = ev.epistemics or {}
                    ev.epistemics['hypothesis'] = True

    # Persist
    store = _load_store()
    store.extend(chosen)
    _save_store(store)

    # Clear pending
    PENDING_PATH.unlink()

    lines = [
        f'Persisted {len(chosen)} statement(s) → statements.json\n'
    ]
    for i, stmt in enumerate(chosen):
        lines.append(_format_stmt_preview(i, stmt))
        lines.append('')

    lines.append(f'Store now contains {len(store)} statements.')
    return '\n'.join(lines)


# ════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', type=int, default=None,
                        help='Run SSE server on this port (omit for stdio)')
    args = parser.parse_args()

    if args.port:
        mcp.run(transport='sse', port=args.port)
    else:
        mcp.run(transport='stdio')
