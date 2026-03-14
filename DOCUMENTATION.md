# Interaction Store — Project Documentation

A persistent, version-controlled framework for curating and querying
mechanistic gene interactions during exploratory research, built on the
[INDRA](https://indra.readthedocs.io) knowledge representation system.

**Status:** Active development — Munch Group, Aarhus University  
**Current focus:** Sex-chromosome gene biology, microtubule-based transport,
archaic introgression analysis

---

## Table of contents

1. [Motivation and design philosophy](#1-motivation-and-design-philosophy)
2. [Repository layout](#2-repository-layout)
3. [Installation and setup](#3-installation-and-setup)
4. [Core concepts](#4-core-concepts)
5. [The statement store](#5-the-statement-store)
6. [INDRA statement types — complete reference](#6-indra-statement-types--complete-reference)
7. [Agent construction](#7-agent-construction)
8. [Evidence construction](#8-evidence-construction)
9. [The gene registry](#9-the-gene-registry)
10. [The MCP server](#10-the-mcp-server)
11. [Claude Code slash commands](#11-claude-code-slash-commands)
12. [Jupyter notebook interface](#12-jupyter-notebook-interface)
13. [Network visualisation](#13-network-visualisation)
14. [INDRA DB enrichment](#14-indra-db-enrichment)
15. [Workflow patterns](#15-workflow-patterns)
16. [Current network — seed content](#16-current-network--seed-content)
17. [Conventions and controlled vocabulary](#17-conventions-and-controlled-vocabulary)

---

## 1. Motivation and design philosophy

Exploratory mechanistic research generates reasoning about gene interactions
long before that reasoning becomes a formal hypothesis, a figure, or a paper.
This framework provides a place to accumulate that reasoning in a structured,
queryable, and version-controlled form from the earliest stages of analysis.

**Three design principles:**

**Semantic precision over free text.**
Every interaction is encoded as a typed INDRA `Statement` rather than a
prose note. This means interactions are queryable by type, by gene, and by
evidence quality; they can be assembled into network graphs; and they can be
cross-referenced against the public INDRA literature database.

**Evidence is first-class.**
Each statement carries one or more `Evidence` objects with provenance —
free-text reasoning, DOI/PMID references, directness annotations, and
hypothesis flags. A statement grows from speculative to literature-supported
by accumulating evidence objects, not by being replaced.

**Separation of concerns.**
Interaction logic (INDRA statements) is kept separate from gene-level
metadata (the registry). The registry carries information that does not
belong in a statement: chromosome, gene group membership, analysis
provenance, expression context, and haplogroup effects. This allows the
same gene to appear in multiple interaction contexts while its metadata
is maintained in one place.

---

## 2. Repository layout

```
project/
│
├── CLAUDE.md                    Reference guide (read by Claude Code on start)
├── .mcp.json                    MCP server config for Claude Code
├── pixi.toml                    Environment and task definitions
│
├── mcp_server.py                FastMCP server — 12 tools
├── gene_registry.py             Registry helpers for notebook use
│
├── statements.json              INDRA statement store  ← auto-managed, commit this
├── gene_registry.json           Gene attribute registry ← auto-managed, commit this
├── pending_extraction.json      Transient — NL extraction staging (do not commit)
│
├── interaction_store.ipynb      Exploratory analysis notebook
├── interaction_network.png      Rendered graph (regenerated on demand)
│
└── .claude/
    └── commands/
        ├── add-interaction.md   /is-add-interaction
        ├── enrich-gene.md       /is-enrich-gene
        ├── rescue-analysis.md   /is-rescue-analysis
        ├── render-network.md    /is-render-network
        └── new-gene-list.md     /is-new-gene-list
```

The two JSON files (`statements.json`, `gene_registry.json`) are the
persistent artefacts. Everything else is tooling around them.
**Both should be committed to version control** — they are the growing
record of the network.

---

## 3. Installation and setup

### Prerequisites

- [Pixi](https://pixi.sh) for environment management
- [Claude Code](https://claude.ai/code) for the agentic interface
- Internet access for INDRA DB queries (optional but recommended)

### Install

```bash
cd project/
pixi install
```

This installs: Python ≥ 3.11, INDRA, FastMCP, NetworkX, Matplotlib,
JupyterLab.

### Start Claude Code with the MCP server

```bash
claude
```

Claude Code reads `.mcp.json` automatically and starts `mcp_server.py`
as a subprocess using stdio transport. The server is available as a
tool set named `interaction-store`.

### Available Pixi tasks

```bash
pixi run mcp          # Start MCP server (stdio — same as Claude Code uses)
pixi run mcp-debug    # Start MCP server on SSE port 8765 for manual testing
pixi run notebook     # Open JupyterLab with the analysis notebook
pixi run render       # Regenerate interaction_network.png from current store
```

### Manual MCP server start (without Pixi)

```bash
python mcp_server.py             # stdio
python mcp_server.py --port 8765 # SSE (debug)
```

---

## 4. Core concepts

### Statement

An INDRA `Statement` is a typed, directed assertion about a relationship
between two biological agents. For example:

```
Phosphorylation(PKA, MAPT, residue='S', position='214')
```

means: PKA phosphorylates MAPT at serine 214. The statement type encodes
the mechanism; the agents encode the participants; optional `residue` and
`position` encode biochemical precision.

### Agent

An `Agent` is a named biological entity (gene, protein, metabolite,
complex) with optional database groundings (`HGNC`, `UniProt`/`UP`,
`CHEBI`, etc.) and optional state conditions (modification state,
mutation, activity).

### Evidence

An `Evidence` object carries the provenance of a statement: the source
(manual reasoning, literature, database), free text, structured
publication identifiers (PMID, DOI), context annotations, and epistemic
flags (hypothesis, negated).

A statement can carry multiple evidence objects. Each additional evidence
object increases the statement's belief score and preserves the
provenance trail from initial hypothesis to confirmed interaction.

### Statement store

`statements.json` — a JSON array of serialised INDRA Statement objects.
Loaded and saved by the store helpers in both `mcp_server.py` and
`interaction_store.ipynb`.

### Gene registry

`gene_registry.json` — a JSON object mapping gene name → attribute dict.
Carries gene-level metadata that does not fit into INDRA's data model:
chromosome, group membership, analysis provenance, rescue logic,
expression context, haplogroup effects, and literature references.

---

## 5. The statement store

### Loading and saving

```python
from indra.statements import stmts_from_json, stmts_to_json
import json, pathlib

def load_store(path='statements.json'):
    p = pathlib.Path(path)
    if not p.exists():
        return []
    with open(p) as f:
        return stmts_from_json(json.load(f))

def save_store(stmts, path='statements.json'):
    with open(path, 'w') as f:
        json.dump(stmts_to_json(stmts), f, indent=2)
```

### Adding a statement

```python
from indra.statements import Phosphorylation, Agent, Evidence

stmt = Phosphorylation(
    Agent('PKA'),
    Agent('MAPT', db_refs={'HGNC': '6893'}),
    residue='S', position='214',
    evidence=[Evidence(
        text='PKA phosphorylates tau at Ser214, reducing MT binding affinity.',
        source_api='manual',
        pmid='12345678',
        text_refs={'PMID': '12345678', 'DOI': '10.1016/j.cell.2020.01.001'},
        annotations={
            'date': '2026-03-11',
            'context': 'cAMP/PKA module',
            'directness': 'direct',
        },
        epistemics={'hypothesis': False},
    )]
)

stmts = load_store()
stmts.append(stmt)
save_store(stmts)
```

### Statement hashing

INDRA computes a shallow hash from statement content (agents, type,
modifications) independent of evidence. This hash is used for
deduplication in preassembly. Adding a second Evidence to an existing
statement is therefore the correct way to promote a hypothesis — not
creating a duplicate statement.

---

## 6. INDRA statement types — complete reference

### 6.1 Activity regulation

Use when the mechanism is unspecified or described only in terms of
functional output.

```python
Activation(subj, obj, obj_activity='activity')
Inhibition(subj, obj, obj_activity='activity')
```

`obj_activity` can specify which activity is regulated:
`'kinase'`, `'phosphatase'`, `'transcription'`, `'gtpbound'`, etc.

**Examples:**
```python
Activation(Agent('IRS2'), Agent('AKT1'))
Inhibition(Agent('ADRA2C'), Agent('PKA'), obj_activity='kinase')
```

### 6.2 Amount regulation

Use when gene expression, mRNA abundance, or protein level is regulated
— not activity.

```python
IncreaseAmount(subj, obj)   # transcriptional activation, stabilisation
DecreaseAmount(subj, obj)   # repression, degradation, silencing
```

**Examples:**
```python
IncreaseAmount(Agent('RBMX2'), Agent('RBMX'))  # paralog backup during MSCI
DecreaseAmount(Agent('PJA1'), Agent('PRKAR1A')) # ubiquitin-mediated degradation
```

> **Activation vs IncreaseAmount:** Use `Activation` for signalling
> cascades where a kinase activates its substrate. Use `IncreaseAmount`
> when the regulated quantity is expression level or abundance.

### 6.3 Post-translational modifications

All PTM types accept optional `residue` (single-letter amino acid) and
`position` (string) arguments. Enzyme is subject; substrate is object.

#### Phosphorylation family

```python
Phosphorylation(enz, sub, residue=None, position=None)
Dephosphorylation(enz, sub, residue=None, position=None)
Autophosphorylation(enz, residue=None, position=None)   # self
Transphosphorylation(enz, residue=None, position=None)  # bound substrate
```

#### Ubiquitin / SUMO / degradation

```python
Ubiquitination(enz, sub, residue=None, position=None)
Deubiquitination(enz, sub, residue=None, position=None)
Sumoylation(enz, sub, residue=None, position=None)
Desumoylation(enz, sub, residue=None, position=None)
```

#### Acetylation / methylation

```python
Acetylation(enz, sub, residue=None, position=None)
Deacetylation(enz, sub, residue=None, position=None)
Methylation(enz, sub, residue=None, position=None)
Demethylation(enz, sub, residue=None, position=None)
```

#### Other covalent marks

```python
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

**Example with site precision:**
```python
Phosphorylation(Agent('PKA'), Agent('MAPT'), residue='S', position='214')
Deacetylation(Agent('HDAC6'), Agent('TUBA1A'), residue='K', position='40')
Ubiquitination(Agent('PJA1'), Agent('PRKAR1A'))
```

### 6.4 Physical association

```python
Complex([agent1, agent2, ...])
```

Use for stable complexes, scaffolding, co-precipitation, and structural
co-localisation where causal direction is absent or irrelevant.

```python
Complex([Agent('AKAP4'), Agent('PKA')])      # scaffolding
Complex([Agent('OFD1'), Agent('SPATC1')])    # co-localisation
```

### 6.5 GTPase regulation

```python
Gef(gef, ras)   # GEF accelerates GDP→GTP exchange (activating)
Gap(gap, ras)   # GAP accelerates GTP hydrolysis (inactivating)
```

### 6.6 Translocation

```python
Translocation(agent, from_location=None, to_location=None)
```

Location strings are GO cellular component terms or common abbreviations:
`'nucleus'`, `'cytoplasm'`, `'plasma membrane'`, `'endosome'`,
`'mitochondria'`, `'axon'`, `'dendrite'`.

```python
Translocation(Agent('NFKB1'), from_location='cytoplasm', to_location='nucleus')
```

### 6.7 Conversion

```python
Conversion(catalyst, inputs=[Agent(...)], outputs=[Agent(...)])
```

Use for metabolic reactions.

```python
Conversion(Agent('ADCY'), inputs=[Agent('ATP')], outputs=[Agent('cAMP')])
```

---

## 7. Agent construction

### Basic

```python
from indra.statements import Agent

Agent('GENE_NAME')
Agent('MAPT', db_refs={'HGNC': '6893', 'UP': 'P10636'})
```

### Convenience wrapper (project standard)

```python
def ag(name: str, hgnc: str = None, up: str = None) -> Agent:
    db_refs = {}
    if hgnc: db_refs['HGNC'] = str(hgnc)
    if up:   db_refs['UP']   = up
    return Agent(name, db_refs=db_refs or None)
```

### With modification state condition

```python
from indra.statements import ModCondition

# Agent that IS phosphorylated at a specific site
Agent('MAPT', mods=[ModCondition('phosphorylation', 'S', '214')])

# Agent that is NOT phosphorylated (negative condition)
Agent('MAPT', mods=[ModCondition('phosphorylation', 'S', '214', is_modified=False)])
```

### With mutation condition

```python
from indra.statements import MutCondition

Agent('BRAF', mutations=[MutCondition('600', 'V', 'E')])  # BRAF V600E
```

### With activity condition

```python
from indra.statements import ActivityCondition

Agent('RAS', activity=ActivityCondition('gtpbound', True))   # active RAS
Agent('GSK3B', activity=ActivityCondition('kinase', False))  # inactive GSK3B
```

### Common HGNC IDs for genes in this network

| Gene | HGNC ID | UniProt |
|---|---|---|
| MAPT | 6893 | P10636 |
| GSK3B | 4617 | P49841 |
| HDAC6 | 14064 | Q9UBN7 |
| MAP7D3 | 28741 | Q8IWZ3 |
| ADRA2C | 283 | P18825 |
| IRS2 | 6136 | O14976 |
| SORCS3 | 17287 | O94933 |
| PRKX | 9441 | P51817 |
| PRKY | 9442 | O43930 |

---

## 8. Evidence construction

### Full Evidence signature

```python
from indra.statements import Evidence

Evidence(
    source_api = 'manual',        # 'manual', 'reach', 'biopax', 'signor', ...
    pmid       = '12345678',      # PubMed ID shortcut (also written to text_refs)
    text       = '...',           # Mechanistic reasoning or supporting sentence
    text_refs  = {                # Structured publication identifiers
        'PMID':  '12345678',
        'DOI':   '10.1016/j.cell.2020.01.001',
        'PMCID': 'PMC7654321',
        'URL':   'https://biorxiv.org/...',
    },
    annotations = {               # Arbitrary key-value metadata
        'date':             '2026-03-11',
        'context':          'cAMP/PKA module',
        'directness':       'direct',       # see §17.1
        'rescue_candidate': True,           # see §17.2
        'species':          'Homo_sapiens', # defaults to Homo_sapiens
    },
    epistemics = {
        'hypothesis': False,      # True if speculative / not yet confirmed
        'negated':    False,      # True if evidence supports ABSENCE of interaction
    },
)
```

### Convenience wrapper (project standard)

```python
from datetime import date

def ev(
    text:       str,
    pmid:       str  = None,
    doi:        str  = None,
    context:    str  = None,
    hypothesis: bool = False,
    direct:     bool = True,
    species:    str  = 'Homo_sapiens',
) -> Evidence:
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
```

### text_refs keys

| Key | Format | Example |
|---|---|---|
| `PMID` | String, no prefix | `'23640493'` |
| `DOI` | Full DOI | `'10.1161/ATVBAHA.113.301522'` |
| `PMCID` | PMC-prefixed | `'PMC3710100'` |
| `URL` | Full URL | `'https://biorxiv.org/content/...'` |

### Multiple evidence objects

Add evidence objects to promote a statement from hypothesis to
literature-supported. The statement hash remains stable; belief score
increases with each supporting evidence.

```python
# Step 1 — initial hypothesis, no paper
stmt = Activation(
    ag('ADRA2C'), ag('PKA'),
    evidence=[ev(
        'ADRA2C suppresses adenylyl cyclase via Gi; net effect is reduced '
        'PKA activity. Inferred from known Gi signalling pathway.',
        context='cAMP/PKA module',
        hypothesis=True,
        direct=False,
    )]
)

# Step 2 — later, paper found; append without replacing
stmt.evidence.append(ev(
    'Haplogroup I males show 0.64-fold lower PRKY expression in macrophages '
    '(P=0.002), consistent with differential cAMP/PKA signalling by haplogroup.',
    pmid='23640493',
    doi='10.1161/ATVBAHA.113.301522',
    context='cAMP/PKA module — Y haplogroup evidence',
    hypothesis=False,
))
save_store(stmts)
```

---

## 9. The gene registry

The registry is a JSON object mapping gene name → attribute dict,
stored in `gene_registry.json`. It is intentionally separate from the
statement store so that gene-level metadata can be queried and updated
independently of the interaction graph.

### Schema

```json
{
  "GENE_NAME": {
    "chromosome": "auto | X | Y | mito",
    "groups": ["group name 1", "group name 2"],
    "analysis_origin": {
      "source":   "IBDmix_NHR | GWAS | eQTL | literature | manual",
      "analysis": "LabName_Method_Year",
      "note":     "Free text: how this gene entered the network"
    },
    "references": [
      {
        "pmid": "12345678",
        "doi":  "10.xxxx/...",
        "note": "Key paper establishing role"
      }
    ],
    "coordinates": {
      "hg38": {"chrom": "chr17", "start": 45894553, "end": 46028334}
    },
    "rescue_logic": "rheostat | paralog_backup | partial | none",
    "contexts": ["neuron", "spermatid", "macrophage"],
    "haplogroup_effect": "Free text describing Y haplogroup-dependent expression",
    "notes": "Free text for anything else"
  }
}
```

### Chromosome values

| Value | Meaning |
|---|---|
| `"auto"` | Autosomal |
| `"X"` | X-linked (subject to MSCI; X-inactivation in females) |
| `"Y"` | Y-linked (haplogroup-variable expression; male-only) |
| `"mito"` | Mitochondrially encoded |

### Python helpers (`gene_registry.py`)

```python
from gene_registry import (
    load_registry, save_registry,
    add_gene, update_gene, add_to_group, add_reference,
    group_members, get_by_chromosome, rescue_candidates,
    enrich_graph, summarise,
)

# Add a new gene
add_gene('PRKY', {
    'chromosome': 'Y',
    'groups': ['cAMP/PKA module'],
    'rescue_logic': 'none',
    'haplogroup_effect': 'Haplogroup I: ~0.64-fold lower expression (Eales 2019)',
    'references': [{'pmid': '23640493', 'doi': '10.1161/ATVBAHA.113.301522'}],
})

# Add a gene to an additional group without overwriting
add_to_group('HDAC6', 'MT lattice/transport')

# Append a reference
add_reference('SORCS3', pmid='29656857', note='SORCS3 in ADHD GWAS')

# Query
x_genes    = get_by_chromosome('X')
mt_members = group_members('MT lattice/transport')
candidates = rescue_candidates()

# Print full summary
summarise()

# Enrich NetworkX graph nodes with registry attributes
enrich_graph(G)   # adds chromosome, groups, rescue_logic to G.nodes
```

### Defined gene groups

| Group | Description |
|---|---|
| `cAMP/PKA module` | ADRA2C, PJA1, AKAP4, PRKX, PRKY, ADCY |
| `IRS2/Akt module` | IRS2, AKT1, GSK3B, TPTE, OCRL, HDAC6 |
| `MT lattice/transport` | MAPT, MAP7D3, HDAC6, DYNLT3 |
| `Neurotrophin/endosome` | SORCS3, DYNLT3, BEX2 |
| `NMD/RNA processing` | XRN2, UPF3B, DDX3X, RBMX2, RBMX, MAF1 |
| `Centriolar/manchette` | SPATC1, OFD1 |
| `NF-κB signalling` | EDA, HIVEP3, NFKB1 |
| `Phosphoinositide` | TPTE, OCRL, IRS2, AKT1 |
| *(analysis-specific groups)* | Created automatically when onboarding gene lists |

---

## 10. The MCP server

`mcp_server.py` is a [FastMCP](https://github.com/jlowin/fastmcp) server
that exposes the statement store and gene registry as tools callable by
Claude Code. It uses stdio transport by default (what Claude Code expects)
and SSE on a configurable port for debugging.

### Configuration (`.mcp.json`)

```json
{
  "mcpServers": {
    "interaction-store": {
      "command": "python",
      "args": ["mcp_server.py"],
      "cwd": "${workspaceFolder}"
    }
  }
}
```

With Pixi:
```json
{
  "mcpServers": {
    "interaction-store": {
      "command": "pixi",
      "args": ["run", "mcp"],
      "cwd": "${workspaceFolder}"
    }
  }
}
```

### Tool reference

#### Statement tools

---

**`add_statement`**

Add a new INDRA statement to the store.

| Parameter | Type | Required | Description |
|---|---|---|---|
| `stmt_type` | str | ✓ | INDRA type name (e.g. `'Phosphorylation'`) |
| `subject` | str | ✓ | Subject gene/protein name |
| `object` | str | ✓ | Object gene/protein name |
| `evidence_text` | str | ✓ | Mechanistic reasoning or supporting text |
| `pmid` | str | | PubMed ID |
| `doi` | str | | DOI |
| `context` | str | | Module tag (e.g. `'cAMP/PKA module'`) |
| `hypothesis` | bool | | `True` if speculative (default `False`) |
| `direct` | bool | | `False` if pathway-level (default `True`) |
| `residue` | str | | Single-letter amino acid for PTMs |
| `position` | str | | Sequence position for PTMs |
| `subject_hgnc` | str | | HGNC ID for subject agent |
| `object_hgnc` | str | | HGNC ID for object agent |

---

**`query_statements`**

Query the store with optional filters.

| Parameter | Type | Description |
|---|---|---|
| `gene` | str | Filter by agent name (subject or object) |
| `stmt_type` | str | Filter by statement type |
| `context` | str | Filter by context annotation (substring match) |
| `hypothesis_only` | bool | Return only hypothetical statements |

---

**`promote_to_literature`**

Append a DOI/PMID-backed `Evidence` to an existing statement, promoting
it from hypothesis to literature-supported.

| Parameter | Type | Description |
|---|---|---|
| `subject` | str | Subject gene name |
| `object` | str | Object gene name |
| `stmt_type` | str | Statement type to match |
| `evidence_text` | str | Supporting text from the paper |
| `pmid` | str | PubMed ID |
| `doi` | str | DOI |
| `context` | str | Context annotation |

---

**`list_statements_summary`**

Returns counts by type, reference coverage, hypothesis count, and all
genes currently in the store.

---

**`enrich_from_indra_db`**

Query the public INDRA DB for known statements connecting a gene list.
Returns a preview — does not auto-persist. Parameters: `genes: list[str]`.

---

#### Registry tools

---

**`add_gene_to_registry`**

Add or update a gene entry with full attribute support.

| Parameter | Type | Description |
|---|---|---|
| `name` | str | Gene symbol |
| `chromosome` | str | `'auto'`, `'X'`, `'Y'`, `'mito'` |
| `groups` | list[str] | Group memberships (appends, does not replace) |
| `analysis_source` | str | `'IBDmix_NHR'`, `'GWAS'`, `'literature'`, etc. |
| `analysis_name` | str | Specific analysis identifier |
| `analysis_note` | str | Free text provenance note |
| `coordinates` | dict | Auto-populated: `{'hg38': {'chrom': ..., 'start': ..., 'end': ...}}` |
| `rescue_logic` | str | `'rheostat'`, `'paralog_backup'`, `'partial'`, `'none'` |
| `contexts` | list[str] | Tissue/cell type list |
| `haplogroup_effect` | str | Y haplogroup-dependent expression note |
| `notes` | str | Free text |
| `pmid` | str | Reference PMID |
| `doi` | str | Reference DOI |
| `reference_note` | str | Note for the reference |

---

**`query_genes`**

Filter the registry.

| Parameter | Type | Description |
|---|---|---|
| `gene` | str | Exact gene name |
| `group` | str | Group membership |
| `chromosome` | str | `'auto'`, `'X'`, `'Y'` |
| `rescue_only` | bool | Return only rescue candidates |
| `analysis` | str | Filter by analysis name |

---

**`get_gene_group`**

Return all members of a named group. Parameter: `group_name: str`.

---

**`list_registry_summary`**

Full registry summary: genes by chromosome, all groups, rescue
candidates, and analysis provenance breakdown.

---

#### NL extraction tools

---

**`extract_from_text`**

Send free text to an NLP reader and preview the extracted statements.
Results are written to `pending_extraction.json` but not persisted.

| Parameter | Type | Description |
|---|---|---|
| `text` | str | One or more mechanistic sentences |
| `reader` | str | `'trips'` (default), `'reach'`, or `'reach_local'` |
| `pmid` | str | PubMed ID to attach to extracted evidence |
| `context` | str | Module/context annotation |
| `service_host` | str | Override default service URL |

Returns a numbered preview. Follow with `persist_extracted` to add chosen
statements.

---

**`persist_extracted`**

Selectively move statements from the pending staging file into the main store.

| Parameter | Type | Description |
|---|---|---|
| `indices` | list[int] | Zero-based indices to persist (all if omitted) |
| `context` | str | Override context annotation |
| `hypothesis` | bool | Mark persisted statements as hypothetical |

Clears `pending_extraction.json` after persisting.

---

#### Browsing tools

---

**`list_contexts`**

List all context tags used in the statement store, with counts.
No parameters. Returns tags sorted by frequency.

---

**`genes_by_context`**

List all genes appearing in statements with a given context tag.

| Parameter | Type | Description |
|---|---|---|
| `context` | str | Context tag (case-insensitive substring match) |

---

**`gene_interactions`**

List all interaction partners for a specific gene, with statement type,
context, and references.

| Parameter | Type | Description |
|---|---|---|
| `gene` | str | Gene symbol |

---

**`gene_info`**

Show everything recorded about a gene: full registry entry (chromosome,
groups, analysis origin, references, notes) plus all interactions from
the statement store.

| Parameter | Type | Description |
|---|---|---|
| `gene` | str | Gene symbol |

---

**`list_groups`**

List all gene groups in the registry with their members. No parameters.

---

#### Graph tools

---

**`render_network`**

Rebuild the interaction graph from the current store and save as PNG.

| Parameter | Type | Description |
|---|---|---|
| `highlight_gene` | str | Dim all nodes except this gene and its direct neighbours |
| `output_path` | str | Output file path (default: `interaction_network.png`) |

---

## 11. Claude Code slash commands

Slash commands live in `.claude/commands/` and are invoked as
`/is-command-name ARGUMENTS`. All interaction store commands use the
`is-` prefix so they are easy to find via tab completion. Each command
is a Markdown prompt template that instructs Claude Code to use the MCP
tools in a specific workflow.

### `/is-add-interaction GENE_A GENE_B`

**Purpose:** Reason about and encode an interaction between two genes.

**Workflow:** Looks up both genes in the registry → selects the most
precise INDRA statement type → assesses directness and rescue candidacy
→ writes evidence text → calls `add_statement` → updates registry if
either gene is new → reports the result.

**Use for:** Any newly identified or reasoned interaction, whether from
literature, database enrichment, or exploratory reasoning.

```
/is-add-interaction IRS2 MAP7D3
/is-add-interaction ADRA2C PJA1
```

---

### `/is-enrich-gene GENE`

**Purpose:** Query INDRA DB for a gene and identify gaps between
database knowledge and the current store.

**Workflow:** Checks current store and registry state → queries INDRA DB
→ summarises top interactions by type and evidence count → identifies
which are not yet in the store → proposes `add_statement` calls →
asks for confirmation before persisting → updates registry.

```
/is-enrich-gene SORCS3
/is-enrich-gene HDAC6
```

---

### `/is-rescue-analysis GENE_A GENE_B`

**Purpose:** Formally assess whether altered expression of one gene
could rescue loss of the other.

**Workflow:** Retrieves current statements and registry data for both
genes → traces the causal chain to any shared downstream effector →
applies the rescue framework (rheostat / paralog / partial / none) →
states the verdict with mechanism and caveats → updates annotations
or appends supporting evidence → compares to existing rescue candidates.

```
/is-rescue-analysis IRS2 MAP7D3
/is-rescue-analysis RBMX2 RBMX
```

---

### `/is-render-network [GENE]`

**Purpose:** Rebuild and save the interaction network graph.

**Workflow:** Reports current store size → calls `render_network` →
reports output path, node/edge counts, and module coverage → identifies
isolated nodes and suggests priority additions.

```
/is-render-network
/is-render-network ADRA2C
```

---

### `/is-extract-from-text TEXT [reader=...] [pmid=...] [context=...]`

**Purpose:** Extract INDRA statements from natural language text, review
them, and selectively persist to the store.

**Workflow:** Sends text to TRIPS or REACH → displays numbered preview of
extracted statements → evaluates each for agent grounding, statement type,
polarity, and directness → recommends which to keep → calls
`persist_extracted` with chosen indices → updates registry for any new
genes → notes any corrections needed for misgrounded agents.

**Reader options:**

| Reader | Endpoint | Best for |
|---|---|---|
| `trips` (default) | trips.ihmc.us | Single sentences, text you write |
| `reach` | agathon.sista.arizona.edu | Sentences from papers |
| `reach_local` | localhost:8080 | Offline; bulk; when remote is down |

```
/is-extract-from-text "ADRA2C inhibits adenylyl cyclase via Gi"
/is-extract-from-text "GSK3β phosphorylates tau at Thr231" reader=reach pmid=12345678
```

---

### `/is-new-gene-list GENE1 GENE2 ...`

**Purpose:** Onboard a batch of genes from a new analysis run.

**Workflow:** Establishes provenance (analysis name, source type, shared
note) → checks existing coverage → classifies each new gene by chromosome
and likely groups → flags SD/pericentromeric concerns as potential
false-positive introgression signals → registers all genes → identifies
priority interactions worth encoding next.

```
/is-new-gene-list HIVEP3 FGGY RASGRP3 HCN1 BTBD3 DOCK2
```

---

### `/is-browse [subcommand] [argument]`

**Purpose:** Browse and list information from the store.

**Subcommands:**

| Subcommand | Description |
|---|---|
| `contexts` | List all context tags with counts |
| `context <TAG>` | List genes in a context |
| `gene <GENE>` | Show full info about a gene (registry + interactions) |
| `interactions <GENE>` | List a gene's interaction partners |
| `groups` | List all gene groups with members |

```
/is-browse contexts
/is-browse context cAMP/PKA
/is-browse gene MAPT
/is-browse interactions PKA
/is-browse groups
```

---

## 12. Jupyter notebook interface

`interaction_store.ipynb` provides an alternative interface to the same
persistent files, suited to longer analytical sessions.

### Cell structure

| Cell | Purpose |
|---|---|
| 1 — Imports | All dependencies; install line at top for Pixi/pip |
| 2 — Helpers | `load_store`, `save_store`, `add_statements`, `ev()`, `ag()`, `summarise()` |
| 3 — Seed | Initial 22-statement network; guard prevents re-seeding |
| 4 — Add new | Template cell for appending new interactions |
| 5 — INDRA DB | Query and optionally persist from public database |
| 6 — Assemble | Build NetworkX graph via `IndraNetAssembler` |
| 7 — Visualise | Kamada-Kawai layout, colour-coded by chromosome |
| 8 — Query | `neighbours()`, `path_between()`, `stmts_for()` helpers |

### Query helpers

```python
# All direct interaction partners of a gene
neighbours('ADRA2C')

# Shortest path between two genes
path_between('ADRA2C', 'MT_lattice')

# All statements involving a gene
stmts = stmts_for('IRS2')
```

### Browsing helpers (`gene_registry.py`)

These functions return Python data structures, suitable for notebook
use, scripting, and programmatic access.

```python
from gene_registry import (
    gene_info,        # full registry entry for one gene
    all_groups,       # all groups with their members
    all_contexts,     # context tags with counts
    genes_by_context, # genes in statements with a context tag → GeneList
    genes_by_group,   # genes belonging to a named group → GeneList
    interactors,      # interaction partners for a gene
    query_statements,     # filter statements by gene/type/context → list[dict]
    query_genes,       # filter registry by gene/group/chrom → dict[str, dict]
)
```

**`gene_info(gene) → dict | None`**

Returns the full registry entry for a gene, or `None` if not registered.

```python
>>> gene_info('MAPT')
{'chromosome': 'auto', 'groups': ['MT lattice/transport'],
 'notes': 'Microtubule-associated protein tau. ...'}
```

**`all_groups() → dict[str, list[str]]`**

Returns `{group_name: [gene, ...]}` for every group in the registry.

```python
>>> all_groups()
{'cAMP/PKA module': ['ADRA2C', 'AKAP4', 'PJA1', 'PRKX', 'PRKY'],
 'gametologs': ['PRKX', 'PRKY'], ...}
```

**`all_contexts(store_path='statements.json') → dict[str, int]`**

Returns `{context_tag: count}` from the statement store.

```python
>>> all_contexts()
{'cAMP/PKA module': 3, 'exploratory': 21, ...}
```

**`genes_by_context(context, store_path='statements.json') → GeneList`**

Returns genes appearing in statements with the given context tag
(case-insensitive substring match) as a `geneinfo.genelist.GeneList`.

```python
>>> genes_by_context('cAMP/PKA')
GeneList(['ADCY', 'ADRA2C', 'AKAP4', 'PKA', 'PJA1', 'PRKX', 'cAMP'])
```

**`genes_by_group(group_name, path=REGISTRY_PATH) → GeneList`**

Returns all genes belonging to a named group as a `geneinfo.genelist.GeneList`.

```python
>>> genes_by_group('MT lattice/transport')
GeneList(['DYNLT3', 'HDAC6', 'MAP7D3', 'MAPT'])
```

**`interactors(gene, store_path='statements.json') → list[dict]`**

Returns all genes that interact with the given gene, with statement type,
context, and references.

```python
>>> interactors('PKA')
[{'gene': 'cAMP', 'type': 'Activation', 'context': 'cAMP/PKA module', 'refs': ['PMID:7803765']},
 {'gene': 'MAPT', 'type': 'Phosphorylation', 'context': 'cAMP/PKA → tau → MT lattice', 'refs': ['PMID:38492709']},
 ...]
```

**`query_statements(intersection=True, **kwargs) → list[dict]`**

Filter the statement store using path-based regex kwargs. Kwarg names
are underscore-separated paths into the nested statement dict — all
possible splits are tried to handle keys that themselves contain
underscores. Values are regexes (`re.search`, case-insensitive).
Multiple filters combine with AND by default; set `intersection=False`
for OR.

```python
>>> query_statements(type='Phospho.*')
[{'type': 'Phosphorylation', 'enz': {'name': 'PKA', ...}, ...}, ...]

>>> query_statements(evidence_text='kinase')         # evidence[*].text
>>> query_statements(subj_name='ADRA2C')             # subj.name
>>> query_statements(evidence_annotations_context='cAMP', type='Activation')
>>> query_statements(evidence_text_refs_PMID='10336') # evidence[*].text_refs.PMID

>>> query_statements(hypothesis_only=True)
[...]  # statements where ≥1 evidence is flagged speculative
```

**`query_genes(intersection=True, **kwargs) → dict[str, dict]`**

Filter the gene registry using path-based regex kwargs. Same path
resolution as `query_statements`. The special kwarg `gene` matches
against the gene name (dict key) rather than a nested path.

```python
>>> query_genes(chromosome='^X$')
{'PRKX': {...}, 'HDAC6': {...}, 'RBMX': {...}, ...}

>>> query_genes(groups='cAMP', chromosome='^X$')  # AND
{'AKAP4': {...}, 'PJA1': {...}, 'PRKX': {...}}

>>> query_genes(analysis_origin_source='IBDmix')
{'ADRA2C': {...}, 'IRS2': {...}, ...}

>>> query_genes(rescue_logic='rheostat')
{'ADRA2C': {...}, 'PJA1': {...}}

>>> query_genes(notes='rescue')
{'ADRA2C': {...}, ...}
```

---

## 13. Network visualisation

### Node colour scheme

| Colour | Meaning |
|---|---|
| Steel blue `#5B9BD5` | Autosomal gene |
| Coral `#D95F5F` | X-linked gene |
| Amber `#E8A020` | Y-linked gene |
| Light grey `#DADDE1` | Signalling intermediate / metabolite |

### Edge colour scheme

| Colour | Statement types |
|---|---|
| Green `#27AE60` | Activation, IncreaseAmount |
| Red `#C0392B` | Inhibition, DecreaseAmount |
| Blue `#2980B9` | Phosphorylation |
| Grey `#888888` | Complex |

### Layout

Kamada-Kawai layout (`nx.kamada_kawai_layout`) with `weight=None` to
treat all edges as equal length. This tends to cluster functional modules
without forcing a strict hierarchy.

### Highlighting a neighbourhood

```python
render_network(highlight_gene='ADRA2C')
# or via MCP:
# render_network(highlight_gene='ADRA2C', output_path='adra2c_neighbourhood.png')
```

Dims all nodes not directly connected to the highlighted gene to 25%
opacity, making local connectivity easier to inspect.

---

## 14. INDRA DB enrichment

The public INDRA database (`db.indra.bio`) contains millions of
pre-assembled statements extracted from the literature by NLP readers
(REACH, Sparser) and curated databases (SIGNOR, PathwayCommons,
BioGRID, PhosphoSitePlus).

### Basic query

```python
from indra.sources import indra_db_rest

# All statements for a gene
stmts = indra_db_rest.get_statements(agents=['SORCS3'], ev_limit=5)

# Directed: what does MAP2K1 phosphorylate?
stmts = indra_db_rest.get_statements(
    subject='MAP2K1', stmt_type='Phosphorylation'
)
```

### Query language (for complex filters)

```python
from indra.sources.indra_db_rest.api import get_statements_from_query
from indra.sources.indra_db_rest.query import HasAgent, HasType, HasDatabases, HasEvidenceBound

# Two genes, database-sourced only, minimum 3 evidence
q = (HasAgent('ADRA2C') & HasAgent('PKA')
     & HasDatabases()
     & HasEvidenceBound(['>= 3']))
result = get_statements_from_query(q)

# Union query
q = HasAgent('MEK') | HasAgent('MAP2K1')
result = get_statements_from_query(q)

# From a specific paper
from indra.sources.indra_db_rest import get_statements_for_paper
stmts = get_statements_for_paper(ids=[('pmid', '12345678')])
```

### Handling DB results

DB results include NLP-extracted statements that may have grounding
errors, polarity errors, or be out of context. Review before persisting:

```python
from indra.statements import pretty_print_stmts
pretty_print_stmts(result.statements[:10])
```

The MCP `enrich_from_indra_db` tool deliberately returns a preview
without auto-persisting for this reason.

---

## 15. Workflow patterns

### Pattern A — Exploratory reasoning → statement

Used during a conversation or analytical session when a new interaction
is identified by reasoning:

1. Identify the mechanism in plain language
2. Choose the most precise statement type (§6)
3. Assess directness (§17.1) and rescue candidacy (§17.2)
4. Write evidence text: name the pathway steps, molecular property
   affected, and any caveats
5. Flag `hypothesis=True` if not yet confirmed by a paper
6. Add via `add_statement` or `add_statements()`
7. Add gene to registry if new

### Pattern B — Literature → statement

Used when a paper directly supports an interaction:

1. Encode the interaction with `hypothesis=False`
2. Include `pmid` and `doi` in the evidence
3. Evidence text should paraphrase the key finding, not quote directly
4. Add modification site (`residue`, `position`) if reported in the paper

### Pattern C — Hypothesis promotion

Used when a paper is found that supports a previously speculative statement:

```python
stmts = load_store()
for s in stmts:
    if (isinstance(s, Activation)
            and s.subj.name == 'ADRA2C'
            and s.obj.name == 'PKA'):
        s.evidence.append(ev(
            'Reduced PRKY expression in haplogroup I men is consistent '
            'with haplogroup-dependent PKA activity differences.',
            pmid='23640493',
            doi='10.1161/ATVBAHA.113.301522',
            hypothesis=False,
        ))
save_store(stmts)
```

### Pattern D — Natural language → statement

Used when you have a precise mechanistic sentence from a paper or from
your own reasoning, and want to let the NLP reader do the initial parsing:

1. Use canonical gene symbols in the text (e.g. GSK3B not GSK-3β)
2. One mechanism per sentence for best accuracy
3. Run extraction: `/is-extract-from-text "sentence" reader=trips`
4. Review the numbered preview — check agent names, type, polarity
5. Persist only the statements that look correct
6. Manually correct any misgrounded agents with `add_statement`

```
/is-extract-from-text "PJA1 ubiquitinates PRKAR2A, targeting it for degradation" pmid=12345678 context="cAMP/PKA module"
```

The two-step `extract_from_text` → `persist_extracted` design is intentional:
NL readers make grounding and polarity errors frequently enough that
auto-persist without review would contaminate the store.

### Pattern E — New analysis run → gene batch

Used when a computational pipeline produces a new list of candidate genes
(e.g. a new IBDmix run, GWAS, eQTL analysis):

```bash
/is-new-gene-list GENE1 GENE2 GENE3 ...
```

The command will ask for provenance details then register all genes with
the analysis name as a group tag, flag any SD/pericentromeric concerns,
and identify priority interactions worth encoding next.

### Pattern F — Module completion

When all major interactions within a module are encoded, run enrichment
to check for known interactions not yet captured:

```bash
/is-enrich-gene SORCS3
/is-enrich-gene DYNLT3
/is-enrich-gene BEX2
```

Then render the network to verify the module is connected:

```bash
/is-render-network SORCS3
```

---

## 16. Current network — seed content

The network was seeded from an exploratory session in March 2026
investigating autosomal–X/Y-linked gene interactions in the context of
microtubule-based transport regulation and MSCI compensation.

### Functional modules

#### cAMP / PKA module

Central axis: `ADRA2C` suppresses adenylyl cyclase (`ADCY`) → ↓cAMP →
reduced activation of `PKA`, `PRKX` (X), and `PRKY` (Y). `PJA1` (X)
opposes this by ubiquitinating PKA regulatory subunits (`PRKAR1A/2A`),
depleting the inhibitory pool and amplifying catalytic PKA activity.
`AKAP4` (X) scaffolds PKA to the fibrous sheath in spermatids.
PKA phosphorylates `MAPT` at Ser214/Thr231.

**Best rescue candidate:** ADRA2C ↔ PJA1 — both modulate the same PKA
catalytic subunit pool by independent mechanisms (cAMP suppression vs.
R-subunit degradation), making them rheostat-like inputs amenable to
mutual compensation.

**Y haplogroup note:** `PRKY` expression is ~0.64-fold lower in
haplogroup I males in macrophages (Eales et al. 2019, PMID 23640493).
This is a cis-regulatory effect on the single-copy Yp11.2 locus, not
CNV — distinct from the ampliconic CNV affecting DAZ/RBMY/TSPY.

#### IRS2 / Akt / MT lattice module

`IRS2` → PI3K → `AKT1` → (1) phospho-inhibits `GSK3B` (Ser9),
reducing tau phosphorylation and releasing MT lattice surface;
(2) phospho-inhibits `HDAC6` (X, Ser22), preserving acetylated tubulin
as the preferred kinesin-1/dynein track.

`MAP7D3` (X) competes with tau for overlapping MT lattice binding sites
and recruits KIF5 for anterograde transport. IRS2 → GSK3β → tau axis
indirectly modulates MAP7D3 binding surface (partial rescue candidate).

`TPTE` is a PTEN paralog that opposes PI3K PIP3 production;
`OCRL` (X) acts on the adjacent PI(4,5)P2 pool (not directly rescuable).

#### Neurotrophin / endosome module

`SORCS3` gates TrkB-containing endosomes to late endosomes after BDNF
binding, controlling retrograde signalling endosome identity. `DYNLT3`
(X) is a Tctex-type dynein light chain required for retrograde transport
of these endosomes. `BEX2` (X) modulates the p75NTR/TrkB co-receptor
interface upstream.

Cargo specification (SORCS3) and motor subunit availability (DYNLT3)
are adjacent but non-redundant steps — not rescuable.

#### NMD / RNA processing module

`UPF3B` (X) generates the 3′ NMD cleavage fragment; `XRN2` degrades it.
Sequential steps — XRN2 upregulation cannot rescue UPF3B loss (wrong
directionality). `DDX3X` (X) cooperates with both in mRNA processing.

`RBMX2` is an autosomal retroposed paralog of `RBMX` (X) that provides
functional backup during MSCI. Both regulate MAPT exon 10 splicing
(3R/4R tau ratio). This is a paralog compensation model.

#### Centriolar / manchette module

`OFD1` (X) controls centriole length; `SPATC1` is centriole-associated
during spermiogenesis. Co-occupancy — no biochemical rescue logic.

#### NF-κB module

`EDA` (X) activates NF-κB via EDAR → TRAF → IKK; `HIVEP3` binds κB
enhancers downstream. Losing EDA signal cannot be rescued by HIVEP3
upregulation — positionally upstream.

### Chromosome breakdown of seed genes

| Chromosome | Genes |
|---|---|
| Autosomal | ADRA2C, IRS2, TPTE, XRN2, RBMX2, SORCS3, SPATC1, HIVEP3, MAPT, GSK3B, AKT1, ADCY, NFKB1 |
| X-linked | PRKX, HDAC6, RBMX, MAP7D3, DYNLT3, UPF3B, DDX3X, PJA1, AKAP4, BEX2, OCRL, OFD1, EDA |
| Y-linked | PRKY |

---

## 17. Conventions and controlled vocabulary

### 17.1 Directness levels

Used in `Evidence.annotations['directness']` and in the `direct`
parameter of the `ev()` wrapper.

| Value | Meaning | Example |
|---|---|---|
| `'direct'` | Protein-protein or sequential enzymatic; adjacent steps | PJA1 ubiquitinates PRKAR1A |
| `'indirect'` | Same pathway, non-adjacent steps | IRS2 → ... → MAP7D3 via tau occupancy |
| `'pathway'` | Same broad process, positionally distinct | BEX2 and SORCS3 on neurotrophin signalling |
| `'structural'` | Co-localisation without biochemical interaction | OFD1 and SPATC1 at centriole |

### 17.2 Rescue logic values

Used in `gene_registry.json['rescue_logic']`.

| Value | Meaning |
|---|---|
| `'rheostat'` | Both genes are independently tunable inputs on the same molecular property — rescue is plausible |
| `'paralog_backup'` | One is a functional paralog of the other on a different chromosome — dosage compensation model |
| `'partial'` | Indirect causal path exists (≥ 3 steps) — partial compensation possible |
| `'none'` | Positionally distinct — different substrates, different branches, or directionality prevents substitution |

### 17.3 Analysis source values

Used in `gene_registry.json['analysis_origin']['source']`.

| Value | Meaning |
|---|---|
| `'IBDmix_NHR'` | NHR candidate from IBDmix introgression analysis |
| `'GWAS'` | Genome-wide association study hit |
| `'eQTL'` | Expression QTL analysis |
| `'literature'` | Identified from manual literature review |
| `'manual'` | Added by hand during exploratory reasoning |

### 17.4 Source API values

Used in `Evidence.source_api`.

| Value | Meaning |
|---|---|
| `'manual'` | Manually curated (this project) |
| `'reach'` | REACH NLP reader |
| `'sparser'` | Sparser NLP reader |
| `'biopax'` | BioPAX database (PathwayCommons) |
| `'signor'` | SIGNOR database |
| `'biogrid'` | BioGRID database |
| `'phosphosite'` | PhosphoSitePlus |
| `'bel'` | BEL (Biological Expression Language) |

### 17.5 Evidence text conventions

- Write in full sentences
- Name the specific pathway steps, not just the endpoint
- Name the molecular property being regulated (e.g. "tau MT binding
  affinity", "PKA R-subunit abundance")
- Include caveats: tissue specificity, directionality constraints,
  cell type restrictions
- For hypothetical statements, make the inferential chain explicit:
  "Inferred from X because Y"
- Maximum ~3 sentences per evidence text

### 17.6 False-positive introgression signal flags

Genes with pericentromeric segmental duplications or known human-specific
SD-embedded copies should be annotated in the registry `notes` field
with the phrase **"potential false-positive introgression signal"** and
a note on the structural reason. Current flagged genes:

- **TPTE** — pericentromeric copies on chr21, chr13, chr15, chr22, Y
- **POTED** — pericentromeric on multiple chromosomes
- **OTOP1 / TMEM128** — pseudocopies in pericentromeric chr2 via
  human-specific SD from 4p16.3
