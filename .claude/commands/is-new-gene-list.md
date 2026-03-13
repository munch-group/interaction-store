# Onboard a new gene list

Onboard the following genes from a new analysis run: $ARGUMENTS

This command is for adding a batch of genes that share a common
analysis provenance (e.g. a new IBDmix run, GWAS hit list, etc.).

Steps:

1. **Establish provenance** — ask the user ONE question at a time.
   Do NOT suggest or guess values — let the user type them.
   After each answer, move to the next question.
   a. "What is the **analysis name**? (format: `LabName_Method_Year`)"
   b. "What is the **source type**? (one of: `IBDmix_NHR`, `GWAS`, `eQTL`, `literature`, `manual`)"
   c. "Any **shared note** or reference (DOI/PMID) for this batch? (or skip)"

   If the user provides all provenance fields in a single message
   (e.g. alongside the gene list or as comma-separated values),
   extract them and confirm — no need to ask again.

2. **Check existing coverage** — call `query_registry` for each gene.
   Report which are already in the registry and which are new.

3. **Classify each new gene**:
   - Chromosome: auto / X / Y
   - Likely gene group(s) based on known biology
   - Whether it has an obvious connection to any gene already in the
     statement store (call `enrich_from_indra_db` in batch if helpful)

4. **Flag structural/SD concerns** — check for genes that are
   known pericentromeric or segmentally duplicated (TPTE, POTED,
   OTOP1/TMEM128 pattern). Flag these as potential false-positive
   introgression signals in the `notes` field.

5. **Register all genes** — call `add_gene_to_registry` for each new
   gene with:
   - `chromosome`
   - `gene_groups` (at minimum `[analysis_name]` as a group)
   - `analysis_source`, `analysis_name`, `analysis_note`

6. **Identify priority interactions** — for genes with clear connections
   to existing store members, list the top 3 interactions worth encoding
   next. Suggest running `/project:add-interaction` for these.

7. **Summary report** — total genes added, chromosome breakdown,
   new groups created, and priority next steps.
