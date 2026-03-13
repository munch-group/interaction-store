# Rescue analysis

Assess whether altered expression of one gene could rescue loss of
the other: $ARGUMENTS

Follow this process:

1. **Retrieve current data** — call `query_statements` and
   `query_registry` for both genes.

2. **Map the causal chain** — trace the mechanistic path from each gene
   to any shared downstream effector. State:
   - The molecular property each gene regulates
   - Whether they regulate the same property or adjacent steps
   - The number of steps between them

3. **Apply the rescue framework** from CLAUDE.md §7:
   - `rheostat`: both are independently tunable inputs on the same
     molecular pool → rescue is plausible
   - `paralog_backup`: one is a functional paralog of the other
     (same biochemical activity, different chromosome) → rescue via
     dosage compensation
   - `partial`: indirect path exists but requires ≥3 steps → partial
     compensation possible
   - `none`: positionally distinct — different substrates, different
     pathway branches, or directionality prevents substitution

4. **State the verdict explicitly**:
   - Rescue plausible? Yes / Partial / No
   - Mechanism: what specifically is being compensated
   - Caveats: tissue specificity, directionality constraints,
     expression context mismatches

5. **Update annotations** — call `promote_to_literature` if there is a
   paper supporting the rescue logic, or update the statement's
   `rescue_candidate` annotation if reasoning alone supports it.

6. **Compare to existing rescue candidates** — call
   `query_registry(rescue_only=True)` and note where this pair ranks
   relative to the current best candidate (ADRA2C ↔ PJA1).
