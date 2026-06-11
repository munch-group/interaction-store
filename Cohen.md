# Verification of X-linked / autosomal mitochondrial gene interactions

This report assesses the 13 STRING-derived interactions between autosomal mitochondrial/RNA-processing genes (MRPS7, GRSF1, PNPT1, PTCD1) and X-linked genes. Each interaction was checked against (1) the INDRA database for high-confidence curated statements and (2) reference verification for the cited literature.

---

## INDRA DB results

Only 3 of 13 interactions have INDRA DB support. All are physical Complex statements with belief scores >= 0.80.

| # | Pair | INDRA type | Belief | Evidence | PMIDs |
|---|------|-----------|--------|----------|-------|
| 2 | GRSF1 -- RPS6KA3 | Complex | 0.82 | 3 (BioGRID) | 32707033, 31678930, 35271311 |
| 3 | GRSF1 -- HSD17B10 | Complex | 0.80 | 4 (3 BioGRID + 1 REACH) | 32877691, 23473034, 29395067 |
| 8 | PTCD1 -- RBMX | Complex | 0.82 | 3 | 26186194, 28514442, 33961781 |

The remaining 10 pairs returned zero INDRA statements. All genes individually exist in INDRA DB with other interaction partners, so the null results are genuine absences rather than lookup failures.

---

## Reference verification

The original report cites 10 references. Four are fully verified; six have incorrect bibliographic metadata (mostly wrong journal names). In all cases the underlying papers exist and the biology is approximately correct.

| # | Cited as | Verdict | Correction |
|---|----------|---------|------------|
| 1 | Menezes et al., 2015, *Hum Mol Genet* (MRPS7 / OXPHOS deficiency) | Verified | PMID 25556185 |
| 2 | Jourdain et al., 2013, *Cell Metabolism* (GRSF1 / MRGs) | Verified | DOI 10.1016/j.cmet.2013.02.005 |
| 3 | Antonicka et al., 2013, *Cell Reports* (GRSF1 / mtRNA processing) | Wrong journal | Actually *Cell Metabolism*, DOI 10.1016/j.cmet.2013.02.006 |
| 4 | He et al., 2018, *Int J Mol Sci* (HSD17B10 / Alzheimer's) | Wrong year | Papers are from 2023, not 2018 (PMID 38139430) |
| 5 | Oerum et al., 2022, *Nature Struct Mol Biol* (HSD17B10 / RNase P) | Wrong journal and year | Actually *J Biol Chem*, 2018 (PMID 29880640) |
| 6 | Perks et al., 2018, *Cell Reports* (PTCD1 / cardiomyopathy) | Exists but misleading | Paper (PMID 29617655) is about ribosome assembly; cardiomyopathy phenotype is in a separate conference abstract (Hughes & Perks et al., 2018, *Heart Lung Circ*) |
| 7 | Fleck et al., 2019, *J Neurosci* (PTCD1 / Alzheimer's) | Verified | PMID 30948477 |
| 8 | Castle et al., 2012, *J Cell Biol* (LAS1L / ITS2 processing) | Wrong journal | Actually *Mol Biol Cell*, 2012 (DOI 10.1091/mbc.e11-06-0530) |
| 9 | Li et al., 2024, *Cancer Res Commun* (LAS1L / SUMOylation) | Verified | PMID 39356143 |
| 10 | PMC12403519, 2025, *J Cell Science* (S6kII / mitochondria in *Drosophila*) | Wrong journal | Actually *Disease Models & Mechanisms* (PMID 40827382) |

The pattern of correct biology but fabricated journal/year details is characteristic of LLM-generated citations.

---

## Interaction-level assessment

### INDRA-backed (strongest support)

**3. GRSF1 -- HSD17B10** (STRING 0.612)
Both proteins co-localize in mitochondrial RNA granules and participate in RNase P-mediated tRNA maturation. INDRA Complex statement with 4 evidence items (belief 0.80). This is the most biologically compelling interaction in the dataset with direct mechanistic evidence.

**2. GRSF1 -- RPS6KA3** (STRING 0.629)
INDRA Complex statement with 3 BioGRID evidence items (belief 0.82). The physical interaction is well-supported. The report's narrative linking RSK2 signaling to mitochondrial dynamics via GRSF1 phosphorylation is plausible but speculative -- no direct evidence for this regulatory axis was found.

**8. PTCD1 -- RBMX** (STRING 0.500)
INDRA Complex statement with 3 evidence items (belief 0.82). Both are RNA-binding proteins; the physical interaction is supported. The proposed regulatory mechanism (RBMX splicing control of PTCD1) is speculative.

### No INDRA support but credible biology

**1. MRPS7 -- EIF1AX** (STRING 0.888)
Highest STRING score in the dataset, driven by experimental (0.581) and coexpression (0.656) channels. Zero INDRA statements. The co-translation coordination argument is biologically reasonable but the framing as a direct interaction is unsupported by curated databases.

**13. PNPT1 -- LAS1L** (STRING 0.410)
Only pair with STRING database evidence (0.360). Zero INDRA statements. The parallel ribosome biogenesis coordination argument is sound.

**4. PNPT1 -- DKC1** (STRING 0.579)
Coexpression (0.411) and textmining (0.268) driven. Zero INDRA statements. Indirect link through parallel ribosome production in different compartments.

**5. GRSF1 -- DDX3X** (STRING 0.510)
Primarily textmining (0.490). Zero INDRA statements despite both genes having large individual statement sets (~700 and ~3000 respectively). The mitochondrial connection narrative is speculative.

**6. PNPT1 -- DDX3X** (STRING 0.506)
Mixed experimental (0.230) and coexpression (0.303). Zero INDRA statements. The innate immunity link (PNPT1 RNA release -> DDX3X/MAVS) is an interesting hypothesis but unverified.

**11. MRPS7 -- UBQLN2** (STRING 0.448)
Experimental (0.302) and coexpression (0.200). Zero INDRA statements. The mitochondrial protein quality control argument is plausible but generic.

### Weakest support

**7. MRPS7 -- KCND1** (STRING 0.503)
No INDRA support. The neuronal co-regulation argument is indirect -- there is no established mechanism linking a mitoribosomal protein to a potassium channel beyond shared tissue expression.

**9. MRPS7 -- NSDHL** (STRING 0.485)
No INDRA support. The MAM contact site narrative is speculative and does not imply a specific functional interaction.

**10. MRPS7 -- SLC35A2** (STRING 0.478)
Driven almost entirely by textmining (0.471). No INDRA support. The connection to mitochondrial biology is tenuous.

**12. MRPS7 -- KDM5C** (STRING 0.437)
Primarily coexpression (0.361). No INDRA support. The epigenetic regulation argument is generic -- KDM5C demethylates thousands of promoters, so coexpression with any nuclear-encoded mitochondrial gene is expected.

---

## Summary

| Tier | Interactions | Basis |
|------|-------------|-------|
| Strong | GRSF1--HSD17B10, GRSF1--RPS6KA3, PTCD1--RBMX | INDRA Complex (belief >= 0.80) + STRING experimental evidence |
| Credible | MRPS7--EIF1AX, PNPT1--LAS1L, PNPT1--DKC1 | High STRING scores but no INDRA support; biology is plausible |
| Speculative | GRSF1--DDX3X, PNPT1--DDX3X, MRPS7--UBQLN2 | Mixed STRING evidence; mitochondrial narratives are hypothetical |
| Weak | MRPS7--KCND1, MRPS7--NSDHL, MRPS7--SLC35A2, MRPS7--KDM5C | No INDRA support; connections to mitochondrial biology are indirect or generic |

Six of 10 cited references have incorrect bibliographic details and should be corrected before the report is used further. HEMK1, noted in the original report as the primary gene of interest in the chr3 region, showed no confident interactions with any X-linked gene.
