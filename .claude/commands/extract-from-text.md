# /project:extract-from-text

Extract INDRA statements from natural language text using an NLP reader,
review the results, then selectively persist to the store.

## Usage

```
/project:extract-from-text TEXT [reader=trips|reach|reach_local] [pmid=XXXXXXXX] [context=MODULE]
```

## Examples

```
/project:extract-from-text "ADRA2C inhibits adenylyl cyclase via Gi, reducing cAMP"
/project:extract-from-text "GSK3β phosphorylates tau at Thr231 and Ser396" reader=reach
/project:extract-from-text "PJA1 ubiquitinates PRKAR2A, targeting it for proteasomal degradation" pmid=12345678 context="cAMP/PKA module"
```

## Workflow

You are helping extract and curate mechanistic interactions from free text.

### Step 1 — Extract

Call `extract_from_text` with the provided text, reader, pmid, and context
arguments. If no reader is specified, use `trips` as the default.

Report the full numbered preview returned by the tool.

### Step 2 — Evaluate each statement

For each extracted statement, reason through:

**Agent grounding**
- Are both agent names correct gene/protein symbols?
- If an agent is named with a common synonym (e.g. "tau" → MAPT,
  "p53" → TP53), note the correct HGNC symbol.
- If grounding looks wrong (e.g. "cAMP" resolved as a protein), flag it.

**Statement type**
- Is the INDRA type correct for the mechanism described?
- Common errors: Activation used where Phosphorylation is more precise;
  IncreaseAmount used where Activation is correct.

**Polarity**
- Is the direction (activation vs inhibition) consistent with the text?
- Negations ("does not activate", "fails to phosphorylate") are a common
  failure mode — check carefully.

**Directness**
- Is this a direct molecular interaction or a pathway-level summary?

**Registry coverage**
- Check whether either agent is already in the registry.
- Note chromosome for any gene not yet registered.

### Step 3 — Recommend which to persist

State clearly which indices you recommend keeping, which to discard, and
why for each. Use this format:

```
KEEP   [0] Phosphorylation(GSK3B, MAPT) — correct type, agents verified
KEEP   [2] Inhibition(ADRA2C, ADCY)     — correct, matches known Gi signalling
DISCARD [1] Activation(cAMP, PKA)       — 'cAMP' is a metabolite not a protein;
                                           encode as Conversion or use Agent with
                                           CHEBI grounding instead
REVIEW [3] Activation(PRKY, PKA)        — direction ambiguous in source text;
                                           flag hypothesis=True
```

### Step 4 — Persist

Call `persist_extracted` with the recommended indices.

For any statement flagged as REVIEW, add `hypothesis=True`.

```python
persist_extracted(indices=[0, 2])
persist_extracted(indices=[3], hypothesis=True)
```

### Step 5 — Registry update

For any genes not yet in the registry, call `add_gene_to_registry` with
at minimum: name, chromosome. Add gene_groups if the interaction fits an
existing module.

### Step 6 — Corrections

If any extracted statement had wrong agent names or statement type,
construct the corrected statement manually using `add_statement` and
note the correction:

```
Note: REACH resolved "tau" as an unknown agent. Corrected manually
to MAPT (HGNC:6893).
```

## Reader selection guide

| Reader | Best for | Endpoint |
|---|---|---|
| `trips` | Single sentences you write; precise mechanistic text | trips.ihmc.us |
| `reach` | Sentences lifted from papers; multi-sentence passages | agathon.sista.arizona.edu |
| `reach_local` | Offline use; bulk processing; when remote is down | localhost:8080 |

Both public endpoints are best-effort research services. If a call fails,
report the error message from the tool and suggest trying the other reader
or `reach_local` with Docker.

## Tips for best extraction quality

- Name both agents explicitly — avoid "it", "this kinase", "the receptor"
- One mechanism per sentence works better than compound sentences
- Use canonical gene symbols (GSK3B not GSK-3β, MAPT not tau)
- Spell out the modification: "phosphorylates" not "modifies"
- Include the site if known: "at Thr231" — REACH captures this as position
