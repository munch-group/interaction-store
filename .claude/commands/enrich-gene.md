# Enrich a gene entry

Enrich the registry and statement store for: $ARGUMENTS

Follow this process:

1. **Check current state** — call `query_registry(gene=...)` and
   `query_statements(gene=...)` to see what is already known.

2. **Query INDRA DB** — call `enrich_from_indra_db([gene])` and summarise
   the top interactions found. Focus on:
   - Statement types with high evidence counts
   - Interactions with other genes already in the store
   - Any Phosphorylation statements with residue/position information

3. **Identify gaps** — compare the DB results with what is already in the
   store. Note any well-supported interactions not yet captured.

4. **Recommend additions** — for each new interaction worth capturing,
   specify the exact `add_statement` call you would make, with full
   evidence text and any available PMID or DOI.

5. **Ask for confirmation** — present the list of proposed additions and
   ask which to persist before calling `add_statement`.

6. **Update registry** — ensure the gene has:
   - Correct `chromosome` field
   - Appropriate `gene_groups`
   - `expression_contexts` if known
   - Any literature references via `add_gene_to_registry`

7. **Summarise** — report what was added and what remains hypothetical.
