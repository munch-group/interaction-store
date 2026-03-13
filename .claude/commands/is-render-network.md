# Render interaction network

Rebuild and save the interaction network graph.

Arguments (optional): $ARGUMENTS
- If a gene name is provided, highlight that gene and its direct neighbours.
- If no argument, render the full network.

Steps:

1. Call `list_statements_summary` to report current store size.

2. Call `render_network` with:
   - `highlight_gene` set to the gene name if one was provided
   - `output_path` left as default unless a path was specified

3. Report:
   - Output file path
   - Node and edge counts
   - Which modules are visible (check registry groups)
   - Any isolated nodes (genes with no connections in the current store)

4. If there are isolated nodes, suggest which interactions would be
   most valuable to add next to connect them.
