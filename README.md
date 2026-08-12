[BnM_PD_Levodopa_LAB_v2.3_README.md](https://github.com/user-attachments/files/30996578/BnM_PD_Levodopa_LAB_v2.3_README.md)
# Genomic screening of probiotic candidates for levodopa compatibility

Data and analysis pipeline for the Perspective *"Why levodopa compatibility should be
genome-screened before probiotic strain selection in Parkinson's disease."*

Everything reported in the paper can be reproduced from what is here. Genome sequences
themselves are **not** included — only the accession list, from which anyone can
re-download the same assemblies.

---

## What was done

5,373 RefSeq assemblies across 22 probiotic species were screened against thresholds
fixed before any genome was retrieved.

| | What it screens | Applied to | Result |
|---|---|---|---|
| **Pre-gate** | resistance and virulence determinants, judged by transferability | 25 *Enterococcus* controls | 25/25 failed |
| **Gate 1** | glycoside hydrolase repertoire | all 5,373 | median 6–31 families by species |
| **Gate 2** | *tyrDC* / *tdc* and related loci | all 5,373 | 570 present (10.6%), 63 ambiguous |

The screen was validated on 27 control genomes **before** the panel was analysed.
All 27 returned the expected call.

---

## Headline results

**Gate 2 is concentrated in four species.**

| Species | n | present | Prevalence (Wilson 95% CI) |
|---|---|---|---|
| *Enterococcus faecium* | 200 | 196 | 98.0% (95.0–99.2) |
| *Enterococcus faecalis* | 201 | 194 | 96.5% (93.0–98.3) |
| *Levilactobacillus brevis* | 200 | 165 | 82.5% (76.6–87.1) |
| *Limosilactobacillus reuteri* | 200 | 15 | 7.5% (4.6–12.0) |
| *Lactiplantibacillus plantarum* | 1,521 | 0 | 0% (upper bound 0.25%) |
| Eighteen further species | 4,572 | 0 | 0% (upper bound 0.084%) |

**A zero is a statement about the calling rule, not about absence of the locus.**
Forty-one *Latilactobacillus curvatus* genomes returned ambiguous calls, and 52 of the
63 ambiguous proteins across the panel are annotated in RefSeq as tyrosine
decarboxylase. `per_genome_calls.tsv` gives every call individually.

**Present calls are not homogeneous.** In *E. faecium*, *E. faecalis* and *L. brevis*
the best target hit is 99.4–100% identical to the characterized reference; in
*L. reuteri* all fifteen fall between 62.8% and 63.4% and are annotated only to family
level.

---

## Contents

```
README.md
run_log.txt                        complete execution record — read this first
pipeline/                          numbered scripts, run in order
data/
  assemblies_all.tsv               all RefSeq assemblies retrieved per taxon
  screened_accessions.txt          the 5,373 accessions actually screened
  gate2_reference.faa              reference proteins with group II decoys
  per_genome_calls.tsv             one row per genome: call, variant, GH count
  table_by_species.tsv             per-species aggregation with Wilson intervals
  gate1_gh_families.txt            GH families detected per genome
  pregate_context.tsv              every resistance/virulence hit in the controls,
                                     with coordinates and adjacent mobile elements
  pregate_per_genome.tsv           pre-gate call for each of the 25 controls
  truncation_controls.faa          synthetic controls used to test variant calling
```

`run_log.txt` records tool versions, retrieval dates, thresholds, the proximity window
used for the pre-gate, and every change made to the analysis after it began.

---

## Reproducing the analysis

```bash
bash pipeline/00_setup.sh                       # datasets, DIAMOND, HMMER
python pipeline/01_reference_proteins.py
python pipeline/02_build_panel.py --taxa taxa_list.txt --out panel/
bash pipeline/03b_download_resume.sh data/screened_accessions.txt
python pipeline/04_validate_controls.py         # STOPS if controls fail
bash pipeline/05_search_all.sh                  # Gate 2
bash pipeline/gate1_fast.sh                     # Gate 1
python pipeline/06_call_gates.py
python pipeline/07_summarize.py
bash pipeline/pregate.sh                        # pre-gate, controls only
python pipeline/pregate_context.py              # mobility context
python pipeline/test_truncation.py              # variant-calling controls
```

Step 04 is not optional. Without it a low observed prevalence cannot be distinguished
from a screen that does not work.

---

## Thresholds, fixed in advance

| Call | Criterion |
|---|---|
| **present** | identity ≥ 60%, coverage ≥ 80%, bitscore margin over next family ≥ 10 |
| **ambiguous** | identity 40–60%, coverage ≥ 50% |
| **absent** | below the above, or highest-scoring hit assigned to a decoy family |

Group II decoys — glutamate (*gadB*), lysine (*ldcA*) and ornithine (*odcA*)
decarboxylases — are searched alongside the target and each protein assigned to its
highest-scoring family. Without them, glutamate decarboxylases are misassigned.

Ambiguous calls are never folded into the negative column. Species with zero detections
are reported with a Wilson upper bound rather than as "0%".

Gate 1 is reported as a distribution, not a pass rate: any cut-off would be arbitrary,
and the gate is conditional on co-culture confirmation.

---

## What the variant classification can and cannot do

`test_truncation.py` builds synthetic controls from a reference TyrDC — truncated to
95, 90, 85, 80, 70 and 50% of length, and with 5, 10, 20 and 40% of residues
substituted — and runs them through the calling rule.

| Control | Call |
|---|---|
| intact | present / intact |
| 5% substituted | present / intact |
| 10–40% substituted | present / **point-variant** |
| truncated to 80–95% | present / intact |
| truncated to 70% or less | **ambiguous** |

**A truncated allele is not resolved under these thresholds.** The coverage requirement
for a present call (≥ 80%) cannot be met by a query short enough to be called
truncated. A screen intended to detect truncation needs a length test applied *before*
the coverage gate. The failure mode this screen detects is substitution.

---

## The pre-gate, on the controls

Applied to the 25 *Enterococcus* control genomes. Between 8 and 54 determinants were
detected per genome, of which 3 to 28 lay within 10 kb of a transposase, integrase,
relaxase or plasmid mobilization function.

**All 123 hits to *van* operon genes were mobile-adjacent, against a background of 46%
for all other determinants (P < 1e-36); at a 5 kb window, 90% and 35%.** The *vanA* and
*vanM* clusters were flanked by IS1216 and IS1542 elements, the *vanB* cluster of
*E. faecalis* V583 by a relaxase and a plasmid mobilization protein.

The background is high because these genomes are IS-dense. An intrinsic efflux pump
present in all 25 genomes was called adjacent in 11 of them, which is indistinguishable
from background. Proximity identifies the clearest cases; it does not partition
determinants cleanly.

---

## Two corrections made during the analysis

Both are in `run_log.txt` with timestamps.

**1. tyrP removed from the calling families.** Tyrosine permease is an operon context
marker, not a call target. Amino acid permeases are ubiquitous and 40–50% similar to
TyrP, so genomes with no *tyrDC* were being called ambiguous. Thresholds were not
changed.

**2. Genomes with no homology hits were being dropped from the denominator.** They are
unambiguous negatives and belong in it. Fixing this changed the overall prevalence from
an inflated 18.6% to 10.6%.

Neither correction was made in response to a result. The first is a coding error; the
second a denominator error.

---

## Known limitations

*Lactiplantibacillus plantarum* IR BL0076 — the one isolate of this species reported to
carry *tyrDC* — has no public genome assembly and is absent from the panel, though the
locus itself is deposited as a sequence (GenBank JQ040309).

Detection of the locus is not evidence of activity. Functional confirmation requires
enzyme assay under controlled tyrosine availability.

---

## Citation

If you use this pipeline or these data, please cite the manuscript. Reference protein
accessions are in Supplementary Table S5; the disclosure audit of published Parkinson's
disease trials is Supplementary Table S3 and S4.

## Licence

Code: MIT. Data: CC-BY-4.0.
