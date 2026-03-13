# Add a mechanistic interaction

Add a new INDRA statement for the genes: $ARGUMENTS

Follow this process:

1. **Identify agents** — extract subject and object gene names from the
   arguments. Look up both in the registry (`query_registry`) to confirm
   chromosome and existing group membership.

2. **Choose statement type** — select the most precise type from CLAUDE.md §1:
   - Prefer modification types (Phosphorylation, Ubiquitination, Acetylation…)
     over Activation/Inhibition when the biochemical mechanism is known.
   - Use Activation/Inhibition when only functional output is described.
   - Use Complex for physical association without clear directionality.
   - Include `residue` and `position` if a modification site is known.

3. **Assess directness** — use CLAUDE.md §7 directness levels:
   - `direct`: protein-protein or sequential enzymatic
   - `indirect`: same pathway, non-adjacent steps
   - `pathway`: same broad process, positionally distinct
   - `structural`: co-localisation without biochemical interaction

4. **Assess rescue candidacy** — reason explicitly about whether altered
   expression of one gene could compensate for loss of the other.
   A rescue is plausible only when both act as rheostats on the same
   molecular property, or when one is a functional paralog of the other.

5. **Flag hypothesis status** — set `hypothesis=True` if the interaction
   is inferred by reasoning rather than directly supported by a paper.

6. **Write the evidence text** — one to three sentences capturing the
   mechanistic reasoning in full. Be precise: name the pathway steps,
   the molecular property affected, and any caveats.

7. **Call `add_statement`** with all parameters populated.

8. **Update the registry** — if either gene is new, call
   `add_gene_to_registry`. If they belong to a new module, add the
   group name.

9. **Report** — summarise: statement type chosen, directness, rescue
   assessment, and whether either gene was newly added to the registry.
