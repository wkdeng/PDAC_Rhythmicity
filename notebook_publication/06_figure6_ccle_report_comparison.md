# Comparison with the supplied CCLE genotyping report

Inspected all three pages of `CCLE Genotyping Data.pdf`. Its 40 calls were compared with the pre-existing Figure S16D evidence table: 33 agree directly, six differ, and one adds an uncertainty qualifier. The PDF links to Cellosaurus and literature but does not identify a DepMap release or supply source mutation/CN rows.


**Reconciliation with supplied `CCLE Genotyping Data.pdf` (3 pages).** The supplied report is a narrative reference summary without an identifiable DepMap release or underlying mutation/CN rows. Its proposed binary matrix is not adopted wholesale. Five display calls are revised using independent source evidence, including the subsequent Panc 10.05 literature refinement:

| Cell line / gene | Previous display | Revised display | Evidence and limitation |
| --- | --- | --- | --- |
| Hs 766T / CDKN2A | WT | mut | [NCI MutSpliceDB](https://brb.nci.nih.gov/cgi-bin/splicing/splicing_evidence.cgi?caid=CA373085801): NM_000077.4:c.457+2T>C with GDC/CCLE RNA-seq splice evidence; homozygosity unverified. |
| Panc 03.27 / CDKN2A | nd | del | [ATCC guide, p10 Table14](https://delivery-files.atcc.org/api/public/content/255660-Cell-Lines-By-Gene-Mutation): homozygous c.1_471del471; retain historical copy-number disagreement. |
| Panc 08.13 / CDKN2A | nd | del | Same ATCC table: homozygous c.1_471del471; reference-stock call. |
| Panc 08.13 / TP53 | WT | nd | [DepMap ACH-000417](https://depmap.org/portal/cell_line/ACH-000417) flags TP53_LoF, conflicting with reported WT; exact alteration unverified. |
| Panc 10.05 / CDKN2A | nd | WT | [2016 Panc10.05 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5097013/) explicitly reports p16(INK4A) WT; [Pa16C study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10010283/) reports CDKN2A WT. |

Panc 03.27 SMAD4 remains **del**: ATCC p28 Table4 reports homozygous c.905_1659del755, contrary to the supplied report's no-clear-alteration summary. Panc 10.05 CDKN2A is **WT**, based on explicit p16/CDKN2A wild-type statements in the [2016 Panc10.05 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5097013/) and the [later Pa16C study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10010283/). [Cellosaurus](https://www.cellosaurus.org/CVCL_1639) lists Pa16C as a synonym. SMAD4 remains **nd (conflicting reports)**: the 2016 study describes a Panc10.05 deletion, whereas the later study describes Pa16C as WT. These are reference annotations; synonymous names do not establish that different laboratories used identical stocks. Current DepMap bulk mutation and CN records were inaccessible during this check, so none of these calls is described as a newly validated current-release DepMap genotype. The source table preserves the report's claim, the previous call, the revised call, and independent evidence for all 40 entries. The supplied report's SHA-256 is `cae30c9c971cb6aea7d6cd8cf7784dd43af429cfa22e504304daee855d189d08`.

Validation: the notebook executed without errors, all 26 inline figure exports were regenerated, five display calls differ from the pre-report annotations (one changed in the Panc 10.05 refinement), and the 17 circadian numerical source tables remained byte-identical.
