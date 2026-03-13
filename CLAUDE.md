# Interaction Store — Reference Guide

This document describes the schema, conventions, and API for the INDRA-based
interaction store used for exploratory mechanistic research in the Munch group.

Three persistent data files live alongside this document:

| File | Contents |
|---|---|
| `statements.json` | INDRA Statement objects (interactions) |
| `gene_registry.json` | Gene/agent metadata and group annotations |
| `indra_db_cache.json` | Per-gene cache of INDRA DB REST API results |

---

## 1. INDRA Statement types

All interactions are encoded as typed INDRA `Statement` objects.
Statements carry **Agents** (subject/object), typed **Evidence**, and
optional site/residue information.

### 1.1 Activity regulation

Use when the mechanism is not fully specified or is at the level of
functional output rather than biochemistry.

```python
Activation(subj, obj)          # subj promotes activity of obj
Inhibition(subj, obj)          # subj suppresses activity of obj
```

Both accept an optional `obj_activity` keyword to specify which activity
is regulated (e.g. `obj_activity='kinase'`).

### 1.2 Post-translational modifications

Enzymatic modifications of a substrate. All accept optional `residue`
(single-letter AA) and `position` (string) arguments for site precision.

```python
# Phosphorylation family
Phosphorylation(enz, sub, residue=None, position=None)
Dephosphorylation(enz, sub, residue=None, position=None)
Autophosphorylation(enz, residue=None, position=None)    # self-modification
Transphosphorylation(enz, residue=None, position=None)   # bound substrate

# Ubiquitin family
Ubiquitination(enz, sub, residue=None, position=None)
Deubiquitination(enz, sub, residue=None, position=None)

# Acetylation / methylation / other marks
Acetylation(enz, sub, residue=None, position=None)
Deacetylation(enz, sub, residue=None, position=None)
Methylation(enz, sub, residue=None, position=None)
Demethylation(enz, sub, residue=None, position=None)
Sumoylation(enz, sub, residue=None, position=None)
Desumoylation(enz, sub, residue=None, position=None)
Hydroxylation(enz, sub, residue=None, position=None)
Dehydroxylation(enz, sub, residue=None, position=None)
Ribosylation(enz, sub, residue=None, position=None)
Deribosylation(enz, sub, residue=None, position=None)
Glycosylation(enz, sub, residue=None, position=None)
Deglycosylation(enz, sub, residue=None, position=None)
Farnesylation(enz, sub, residue=None, position=None)
Defarnesylation(enz, sub, residue=None, position=None)
Geranylgeranylation(enz, sub, residue=None, position=None)
Degeranylgeranylation(enz, sub, residue=None, position=None)
Palmitoylation(enz, sub, residue=None, position=None)
Depalmitoylation(enz, sub, residue=None, position=None)
Myristoylation(enz, sub, residue=None, position=None)
Demyristoylation(enz, sub, residue=None, position=None)
```

Example with site precision:
```python
Phosphorylation(ag('PKA'), ag('MAPT'), residue='S', position='214')
```

### 1.3 Amount regulation

Use when expression level / abundance (not activity) is being regulated.

```python
IncreaseAmount(subj, obj)      # subj increases expression/abundance of obj
DecreaseAmount(subj, obj)      # subj decreases expression/abundance of obj
```

**Distinction from Activation/Inhibition**: use `IncreaseAmount`/`DecreaseAmount`
for transcriptional regulation, mRNA stability, protein degradation.
Use `Activation`/`Inhibition` for signalling/activity effects.

### 1.4 Physical association

```python
Complex([agent1, agent2, ...])  # stable complex; no directionality
```

Use for scaffolding, co-precipitation, structural co-localisation where
causal direction is unknown or irrelevant.

### 1.5 GTPase regulation

```python
Gef(gef, ras)   # GEF accelerates GDP→GTP exchange (activating)
Gap(gap, ras)   # GAP accelerates GTP hydrolysis (inactivating)
```

### 1.6 Translocation

```python
Translocation(agent, from_location=None, to_location=None)
```

Use for nuclear/cytoplasmic shuttling, membrane recruitment, etc.
Locations are GO cellular component terms or common abbreviations
(`'nucleus'`, `'cytoplasm'`, `'plasma membrane'`, `'endosome'`).

### 1.7 Conversion

```python
Conversion(catalyst, inputs=[Agent(...)], outputs=[Agent(...)])
```

Use for metabolic/enzymatic reactions that convert substrates to products.

---

## 2. Agent construction

```python
from indra.statements import Agent

# Minimal
Agent('GENE_NAME')

# With database grounding (preferred for INDRA DB queries)
Agent('MAPT', db_refs={'HGNC': '6893', 'UP': 'P10636'})

# With modification state condition
from indra.statements import ModCondition
Agent('MAPT', mods=[ModCondition('phosphorylation', 'S', '214')])

# With mutation condition
from indra.statements import MutCondition
Agent('BRAF', mutations=[MutCondition('600', 'V', 'E')])  # V600E

# With activity condition
from indra.statements import ActivityCondition
Agent('RAS', activity=ActivityCondition('gtpbound', True))
```

Convenience wrapper used in this project:
```python
def ag(name, hgnc=None, up=None):
    db_refs = {}
    if hgnc: db_refs['HGNC'] = str(hgnc)
    if up:   db_refs['UP'] = up
    return Agent(name, db_refs=db_refs or None)
```

---

## 3. Evidence construction

`Evidence` objects carry provenance for each statement.
Multiple `Evidence` objects on one statement mean multiple sources
support the same interaction — adding a second `Evidence` with a
`pmid`/`doi` promotes a manual hypothesis to literature-grounded.

### Full Evidence signature

```python
Evidence(
    source_api = 'manual',         # 'manual', 'reach', 'biopax', 'signor', ...
    pmid       = '12345678',       # PubMed ID shortcut (also goes into text_refs)
    text       = 'mechanism ...',  # supporting sentence or reasoning
    text_refs  = {                 # structured publication identifiers
        'PMID':  '12345678',
        'DOI':   '10.1161/ATVBAHA.113.301522',
        'PMCID': 'PMC1234567',
        'URL':   'https://...',
    },
    annotations = {                # arbitrary key-value metadata
        'date':             '2026-03-11',
        'context':          'cAMP/PKA module',
        'directness':       'direct',   # see §7
        'rescue_candidate': True,
    },
    epistemics = {
        'hypothesis': False,       # True if speculative / not yet confirmed
        'negated':    False,       # True if evidence supports ABSENCE of interaction
    },
)
```

### text_refs keys

| Key | Format | Example |
|---|---|---|
| `PMID` | string | `'23640493'` |
| `DOI` | string | `'10.1161/ATVBAHA.113.301522'` |
| `PMCID` | string | `'PMC3710100'` |
| `URL` | string | `'https://www.biorxiv.org/...'` |

### Convenience wrapper (current)

```python
def ev(text, pmid=None, doi=None, context=None, hypothesis=False, direct=True):
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
            'context':    context,
            'directness': 'direct' if direct else 'indirect',
        },
        epistemics={'hypothesis': hypothesis},
    )
```

### Example: seeding then promoting a hypothesis

```python
# Step 1 — add as hypothesis (no paper yet)
Activation(
    ag('ADRA2C'), ag('PKA'),
    evidence=[ev(
        'ADRA2C suppresses cAMP; net effect reduces PKA activity.',
        context='cAMP/PKA module',
        hypothesis=True,
        direct=False,
    )]
)

# Step 2 — later, promote with literature reference by appending Evidence
stmt.evidence.append(ev(
    'Haplogroup I carriers show reduced UTY and PRKY expression in macrophages.',
    pmid='23640493',
    doi='10.1161/ATVBAHA.113.301522',
    context='cAMP/PKA module — Y haplogroup effect',
    hypothesis=False,
))
```

---

## 4. Gene registry

A sidecar `gene_registry.json` stores gene-level attributes that do not
belong in INDRA statements. The schema is:

```json
{
  "GENE_NAME": {
    "chromosome": "auto | X | Y | mito",
    "gene_groups": ["cAMP/PKA module", "MT transport"],
    "analysis_origin": {
      "source":   "IBDmix_NHR | manual | GWAS | literature",
      "analysis": "Tishkoff_IBDmix_2026 | ...",
      "note":     "Free text describing how this gene entered the network"
    },
    "references": [
      {
        "pmid":  "12345678",
        "doi":   "10.xxxx/...",
        "note":  "Key paper establishing this gene's role"
      }
    ],
    "rescue_logic": "rheostat | paralog_backup | pathway_neighbor | none",
    "expression_contexts": ["neuron", "spermatid", "macrophage"],
    "haplogroup_effect": "haplogroup I: 0.64-fold lower PRKY expression (Eales 2019)",
    "notes": "Free text"
  }
}
```

### Built-in chromosome groups

| Value | Meaning |
|---|---|
| `"auto"` | Autosomal gene |
| `"X"` | X-linked (subject to MSCI, dosage compensation) |
| `"Y"` | Y-linked (haplogroup-variable expression) |
| `"mito"` | Mitochondrially encoded |

### Built-in gene groups (current)

| Group | Members |
|---|---|
| `cAMP/PKA module` | ADRA2C, PJA1, AKAP4, PRKX, PRKY |
| `IRS2/Akt module` | IRS2, AKT1, GSK3B, TPTE, OCRL, HDAC6 |
| `MT lattice/transport` | MAPT, MAP7D3, HDAC6, DYNLT3 |
| `Neurotrophin/endosome` | SORCS3, DYNLT3, BEX2 |
| `NMD/RNA processing` | XRN2, UPF3B, DDX3X, RBMX2, RBMX |
| `Centriolar/manchette` | SPATC1, OFD1 |
| `NF-κB signalling` | EDA, HIVEP3 |
| `Phosphoinositide` | TPTE, OCRL, IRS2 |
| `IBDmix NHR candidates` | HIVEP3, FGGY, RASGRP3, SORCS3, MYO16, IRS2, BTBD3, XRN2, NKX2-4, RBMX2, HCN1, TPTE, POTED, ADRA2C, OTOP1, TMEM128, MAF1, SPATC1 |

---

## 5. Registry helpers (gene_registry.py)

```python
from gene_registry import load_registry, save_registry, add_gene, get_group

# Load
reg = load_registry()           # dict: gene_name → attributes

# Add / update a gene
add_gene('PRKY', {
    'chromosome': 'Y',
    'gene_groups': ['cAMP/PKA module'],
    'analysis_origin': {
        'source': 'literature',
        'note': 'Y-linked PKA paralog; haplogroup I shows reduced expression'
    },
    'references': [{'pmid': '23640493', 'note': 'Eales 2019 - UTY/PRKY in macrophages'}],
    'haplogroup_effect': 'haplogroup I: ~0.64-fold lower expression (P=0.002)',
})

# Query by group
mt_genes = get_group('MT lattice/transport')

# Enrich NetworkX graph with registry attributes
for node in G.nodes():
    info = reg.get(node, {})
    G.nodes[node]['chromosome']  = info.get('chromosome', 'unknown')
    G.nodes[node]['gene_groups'] = info.get('gene_groups', [])
```

---

## 6. INDRA DB queries and caching

### 6.1 Cache layer (`indra_cache.py`)

All INDRA DB REST queries should go through the persistent cache in
`indra_cache.py`. Results are stored per-gene in `indra_db_cache.json`
so repeated queries are instant instead of taking minutes per gene.

```python
from indra_cache import (
    cached_get_statements,        # single gene
    cached_get_statements_batch,  # multiple genes with progress
    cache_summary,                # show cached genes and ages
    invalidate_cache,             # clear specific genes or all
)

# Single gene — checks cache, queries API on miss
stmts = cached_get_statements('ADRA2C', ev_limit=5)

# Batch — prints (cached) or (fetched) per gene
all_stmts = cached_get_statements_batch(
    ['ADRA2C', 'IRS2', 'MAPT'],
    ev_limit=5,
    max_age_days=30,   # re-fetch entries older than 30 days (None = never expire)
)

# Inspect cache state
print(cache_summary())

# Force refresh for specific genes
invalidate_cache(['ADRA2C', 'IRS2'])

# Clear entire cache
invalidate_cache()
```

The MCP `enrich_from_indra_db` tool and the notebook both use this
cache automatically.

### 6.2 Direct INDRA DB REST API

For queries not covered by the cache (directed queries, paper lookups),
use the INDRA API directly. Note that `get_statements()` returns a
processor object — access `.statements` for the list.

```python
from indra.sources.indra_db_rest import get_statements
from indra.sources.indra_db_rest.api import get_statements_from_query
from indra.sources.indra_db_rest.query import HasAgent, HasType, HasDatabases

# All statements for a gene (returns processor, not list)
proc = get_statements(agents=['ADRA2C'], ev_limit=5)
stmts = proc.statements

# Directed query: what does MAP2K1 phosphorylate?
proc = get_statements(subject='MAP2K1', stmt_type='Phosphorylation')
stmts = proc.statements

# Two-gene query with evidence count filter
q = HasAgent('ADRA2C') & HasAgent('PKA') & HasEvidenceBound(['>= 2'])
result = get_statements_from_query(q)

# Database-sourced only (no NLP extraction noise)
q = HasAgent('SORCS3') & HasDatabases()
result = get_statements_from_query(q)

# From a specific paper
from indra.sources.indra_db_rest import get_statements_for_paper
stmts = get_statements_for_paper(ids=[('pmid', '12345678')])
```

---

## 7. Workflow conventions

### Adding an interaction from exploratory reasoning (no paper yet)
```python
Activation(
    ag('GENE_A'), ag('GENE_B'),
    evidence=[ev(
        'Mechanism text capturing the reasoning in full.',
        context='Module name',
        hypothesis=True,    # flag as speculative if not yet confirmed
        direct=False,       # flag if pathway-level, not protein-protein
    )]
)
```

### Promoting a hypothesis to literature-supported
Append a second `Evidence` object with `pmid`/`doi` populated.
Keep the original manual evidence — multiple evidence objects
raise the statement's belief score in INDRA's preassembly and
preserve the reasoning trail.

```python
# Find the statement and append
stmts = load_store()
for s in stmts:
    if isinstance(s, Activation) and \
       s.subj.name == 'GENE_A' and s.obj.name == 'GENE_B':
        s.evidence.append(ev(
            'Direct quote or paraphrase from paper.',
            pmid='12345678',
            doi='10.1000/xyz123',
            context='Module name',
            hypothesis=False,
        ))
save_store(stmts)
```

### Directness levels (use in `annotations['directness']`)

| Value | Meaning |
|---|---|
| `'direct'` | Protein-protein or sequential enzymatic (e.g. PJA1 ubiquitinates PRKAR1A) |
| `'indirect'` | Same pathway, non-adjacent steps (e.g. IRS2 → ... → MAP7D3 via tau) |
| `'pathway'` | Same broad process, positionally distinct (e.g. BEX2 / SORCS3 on neurotrophin) |
| `'structural'` | Co-localisation without biochemical interaction (e.g. OFD1 / SPATC1) |

### Rescue candidacy (use in `annotations['rescue_candidate']`)

| Value | Meaning |
|---|---|
| `True` | Altered expression of one gene could plausibly compensate for loss of the other |
| `'partial'` | Partial compensation possible (longer causal chain) |
| `False` | Positionally distinct — rescue not plausible |

Best current rescue candidate: **ADRA2C ↔ PJA1** (both rheostat-like
inputs on the same PKA R-subunit pool via independent mechanisms).

---

## 8. MCP server tools

The local MCP server (`mcp_server.py`) exposes these tools to Claude Code:

| Tool | Purpose |
|---|---|
| `add_statement` | Add a typed INDRA statement with full evidence |
| `query_statements` | Filter store by gene, type, context, hypothesis flag |
| `promote_to_literature` | Append a DOI/PMID evidence to an existing statement |
| `list_statements_summary` | Counts by type, reference coverage, genes covered |
| `enrich_from_indra_db` | Query public INDRA DB for a gene list (cached, preview only) |
| `add_gene_to_registry` | Add/update gene with chromosome, groups, provenance, refs |
| `query_registry` | Filter registry by gene, group, chromosome, rescue logic |
| `get_gene_group` | List all members of a named group |
| `list_registry_summary` | Full registry: chromosomes, groups, rescue candidates |
| `render_network` | Rebuild and save the interaction graph PNG |

### Starting the server

```bash
# Claude Code reads .mcp.json automatically — server starts on `claude`
pixi run mcp        # stdio (production)
pixi run mcp-debug  # SSE on port 8765 (debugging)
```

---

## 9. Slash commands

| Command | Usage |
|---|---|
| `/project:add-interaction GENE_A GENE_B` | Reason about and encode an interaction |
| `/project:enrich-gene GENE` | Query INDRA DB and update registry |
| `/project:rescue-analysis GENE_A GENE_B` | Assess rescue candidacy |
| `/project:render-network [GENE]` | Rebuild graph, optionally highlighting a gene |
| `/project:new-gene-list GENE1 GENE2 ...` | Onboard a batch from an analysis run |

---

## 10. File layout

```
project/
├── CLAUDE.md                   ← this file (read automatically by Claude Code)
├── .mcp.json                   ← MCP server config for Claude Code
├── pixi.toml                   ← environment and task definitions
├── mcp_server.py               ← FastMCP server (10 tools)
├── statements.json             ← INDRA statement store (auto-managed)
├── gene_registry.json          ← gene attribute registry (auto-managed)
├── indra_db_cache.json         ← per-gene INDRA DB query cache (auto-managed)
├── gene_registry.py            ← registry helpers for notebooks
├── indra_cache.py              ← INDRA DB cache layer (used by MCP + notebook)
├── interaction_store.ipynb     ← exploratory analysis notebook
├── interaction_network.png     ← rendered graph (regenerated on demand)
└── .claude/
    └── commands/
        ├── add-interaction.md
        ├── enrich-gene.md
        ├── rescue-analysis.md
        ├── render-network.md
        └── new-gene-list.md
```
