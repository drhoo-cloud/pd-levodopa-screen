# Genomic screening of lactic acid bacteria for levodopa compatibility

Data and analysis pipeline for the manuscript *"A sequence-first qualification
framework for levodopa-compatible probiotic candidates in Parkinson's disease"*.

Everything reported in the paper can be reproduced from what is in this
repository. Genome sequences themselves are **not** included — only the
accession list, from which anyone can re-download the same assemblies.

---

## What was done

5,373 RefSeq genome assemblies across 22 lactic acid bacterial and
bifidobacterial species were screened against two loci, using thresholds fixed
before any genome was retrieved.

| | What it screens | Result |
|---|---|---|
| **Gate 2** | *tyrDC* / *tdc* and related group II decarboxylase loci | 570 present (10.6%), 63 ambiguous |
| **Gate 1** | Glycoside hydrolase (GH) family repertoire | median 6–31 families per genome by species |

The screen was validated on 27 control genomes **before** the panel was
analysed. All 27 returned the expected call.

---

## Headline results

**Gate 2 is concentrated in four species and absent from the other eighteen.**

| Species | n | *tyrDC*/*tdc* present | Prevalence (Wilson 95% CI) |
|---|---|---|---|
| *Enterococcus faecium* | 200 | 196 | 98.0% (95.0–99.2) |
| *Enterococcus faecalis* | 201 | 194 | 96.5% (93.0–98.3) |
| *Levilactobacillus brevis* | 200 | 165 | 82.5% (76.6–87.1) |
| *Limosilactobacillus reuteri* | 200 | 15 | 7.5% (4.6–12.0) |
| *Lactiplantibacillus plantarum* | 1,521 | 0 | 0% (upper bound 0.25%) |
| Seventeen further species (pooled) | 3,051 | 0 | 0% (upper bound 0.13%) |

**The two gates are independent once genus is held constant.** Within species,
carriers and non-carriers have the same GH repertoire: *E. faecalis* 31 vs 31,
*L. brevis* 14 vs 15, *L. reuteri* 11 vs 11.

---

## Contents

```
README.md                          this file
run_log.txt                        complete execution record — read this first
pipeline/                          numbered scripts, run in order
data/
  assemblies_all.tsv               all RefSeq assemblies retrieved per taxon
  screened_accessions.txt          the 5,373 accessions actually screened
  gate2_reference.faa              44 reference proteins, labelled by family
  per_genome_calls.tsv             one row per genome: Gate 2 call, variant, GH count
  table_by_species.tsv             per-species aggregation with Wilson intervals
  gate1_gh_families.txt            GH families detected per genome
```

`run_log.txt` records tool versions, retrieval dates, thresholds, genome counts,
and every change made to the analysis after it began — including the two
corrections described below.

---

## Reproducing the analysis

```bash
bash pipeline/00_setup.sh                       # install datasets, DIAMOND, HMMER
python pipeline/01_reference_proteins.py        # fetch reference proteins
python pipeline/02_build_panel.py --taxa taxa_list.txt --out panel/
bash pipeline/03b_download_resume.sh data/screened_accessions.txt
python pipeline/04_validate_controls.py         # STOPS if controls fail
bash pipeline/05_search_all.sh                  # Gate 2
bash pipeline/gate1_fast.sh                     # Gate 1
python pipeline/06_call_gates.py
python pipeline/07_summarize.py
```

Step 04 is not optional. Without it a low observed prevalence cannot be
distinguished from a screen that does not work.

---

## Thresholds, fixed in advance

| Call | Criterion |
|---|---|
| **present** | identity ≥ 60%, coverage ≥ 80%, bitscore margin over next family ≥ 10 |
| **ambiguous** | identity 40–60%, coverage ≥ 50% |
| **absent** | below the above, or highest-scoring hit assigned to a decoy family |

Group II decoys — glutamate (*gadB*), lysine (*ldcA*) and ornithine (*odcA*)
decarboxylases — are searched alongside the target and each protein is assigned
to its highest-scoring family. Without them, glutamate decarboxylases are
misassigned to the target set.

Ambiguous calls are never folded into the negative column. Species with zero
detections are reported with a Wilson upper bound rather than as "0%".

Gate 1 is reported as a distribution, not a pass rate. No GH threshold was
imposed: any cut-off would be arbitrary at this stage, and the gate is in any
case conditional on co-culture confirmation.

---

## Two corrections made during the analysis

Both are in `run_log.txt` with timestamps.

**1. tyrP removed from the calling families.** Tyrosine permease is an operon
context marker, not a call target. Amino acid permeases are ubiquitous and
40–50% similar to TyrP, so genomes with no *tyrDC* were being called ambiguous.
Thresholds were not changed.

**2. Genomes with no homology hits were being dropped from the denominator.**
They are unambiguous negatives and belong in it. Fixing this changed the overall
prevalence from an inflated 18.6% to 10.6%.

Neither correction was made in response to a result. The first is a coding
error; the second is a denominator error.

---

## Known limitations

*Lactiplantibacillus plantarum* IR BL0076 — the one isolate of this species
reported to carry *tyrDC*, at 98% identity to that of *L. brevis* — has no public
assembly and is therefore absent from the panel. The 0/1,521 result should be
read with that in mind, and it is itself an instance of the disclosure limit the
paper documents.

Detection of the locus is not evidence of activity. Every positive call is
classified as intact, truncated or point-variant, but functional confirmation
requires enzyme assay under controlled tyrosine availability.

---

## Citation

If you use this pipeline or these data, please cite the manuscript. Reference
protein accessions are listed in Supplementary Table S7; the complete audit of
strain disclosure in published Parkinson's disease trials is Supplementary
Table S6.

## Licence

Code: MIT. Data: CC-BY-4.0.
