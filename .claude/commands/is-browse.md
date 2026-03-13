# Browse the interaction store

Browse and list information from the store. Usage:

- `/is-browse contexts` — list all context tags with counts
- `/is-browse context <TAG>` — list genes in a context
- `/is-browse gene <GENE>` — show all info about a gene
- `/is-browse interactions <GENE>` — list a gene's interactions
- `/is-browse groups` — list all gene groups

Arguments: $ARGUMENTS

---

Parse the arguments and call the appropriate MCP tool:

- "contexts" or no arguments → call `list_contexts`
- "context <TAG>" → call `genes_by_context(context=TAG)`
- "gene <GENE>" or just a gene name → call `gene_info(gene=GENE)`
- "interactions <GENE>" → call `gene_interactions(gene=GENE)`
- "groups" → call `list_groups`

Present results directly — no extra commentary needed.
