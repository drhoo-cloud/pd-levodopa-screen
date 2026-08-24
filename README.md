# Genomic screening of probiotic candidates for levodopa compatibility

Data and analysis pipeline for the Perspective *"Genome-screened probiotic selection for
levodopa compatibility and gut–brain axis support in Parkinson's disease."*

Everything reported in the paper can be reproduced from what is here. Genome sequences
themselves are **not** included — only the accession list, from which anyone can
re-download the same assemblies.

> **Revised 24 August 2026.** Both gates were re-run after peer review of the analysis
> by a co-author. Gate 1 now requires agreement between at least two of the three dbCAN
> prediction methods; Gate 2 now scores each genome against a profile HMM built from
> UniProt-reviewed bacterial tyrosine decarboxylases rather than a DIAMOND search against
> a reference set containing decoys. Numbers below are the revised ones. See
> [Revision history](#revision-history).

---

## What was done

5,373 RefSeq assemblies across 22 probiotic species were screened against thresholds
fixed before any genome was retrieved.

| | What it screens | Applied to | Result |
|---|---|---|---|
| **Pre-gate** | resistance and virulence determinants, judged by transferability; histidine decarboxylase | 25 *Enterococcus* controls | 25/25 failed |
| **Gate 1** | glycoside hydrolase repertoire | all 5,373 | species medians 7–27 families |
| **Gate 2** | *tyrDC* / *tdc* | all 5,373 | 633 present (11.8%) |

The screen was validated on 27 control genomes **before** the panel was analysed.
All 27 returned the expected call.

---

## Headline results

**Gate 2 is concentrated in six species.**

| Species | n | present | Prevalence (Wilson 95% CI) |
|---|---|---|---|
| *Enterococcus faecium* | 200 | 200 | 100.0% (98.1–100.0) |
| *Enterococcus faecalis* | 201 | 194 | 96.5% (93.0–98.3) |
| *Levilactobacillus brevis* | 200 | 166 | 83.0% (77.2–87.6) |
| *Latilactobacillus curvatus* | 117 | 41 | 35.0% (27.0–44.0) |
| *Limosilactobacillus reuteri* | 200 | 28 | 14.0% (9.9–19.5) |
| *Pediococcus pentosaceus* | 200 | 4 | 2.0% (0.8–5.0) |
| *Lactiplantibacillus plantarum* | 1,521 | 0 | 0% (upper bound 0.25%) |
| Sixteen further species | 4,255 | 0 | 0% (upper bound 0.09%) |

**The score distribution is bimodal, which is why there is no ambiguous column.**
Negative genomes score between 0% and 4.3% of the profile self-score; positives between
61.5% and 100%. Exactly one assembly falls in between, at 19.5%. Any cut-off between
roughly 20% and 60% returns the same 633 calls. The earlier version of this screen
produced 63 ambiguous calls from a bitscore-margin rule; all of them resolve cleanly
under the normalized score.

**Present calls are not homogeneous.** Normalized profile scores among positive genomes:

| Species | Score range |
|---|---|
| *Enterococcus faecalis* | 99.8–100.0% |
| *Enterococcus faecium* | 93.0% |
| *Levilactobacillus brevis* | 68.4–96.1% |
| *Latilactobacillus curvatus* | 61.5–88.3% |
| *Limosilactobacillus reuteri* | 71.4–77.1% |
| *Pediococcus pentosaceus* | 61.5–62.0% |

Three of the four reference sequences are enterococcal, so a low-but-positive score
reflects distance from the reference set as well as divergence at the locus. Calls at
the bottom of the positive mode are referrals to enzyme assay, not graded activity.
`per_genome_calls.tsv` gives every call individually.

---

## Contents

```
README.md
run_log.txt                        complete execution record — read this first
APPLY.md                           what changed in the August 2026 revision, and why
pipeline/                          numbered scripts, run in order
refs/
  tyrdc_reference.faa              the four UniProt-reviewed reference sequences
  tyrdc.aln                        MAFFT alignment
  tyrdc.hmm                        profile HMM built from the alignment
  self_score.txt                   1342.2 — the denominator for every normalized score
  gate2_profile.json               query, release, model length, MD5, self-score
  pregate_hdc.faa                  histidine decarboxylase references for the pre-gate
data/
  assemblies_all.tsv               all RefSeq assemblies retrieved per taxon
  screened_accessions.txt          the 5,373 accessions actually screened
  per_genome_calls.tsv             one row per genome: call, score, variant, GH count
  table_by_species.tsv             per-species aggregation with Wilson intervals
  gate1_gh_families.txt            GH families detected per genome (2-of-3 consensus)
  pregate_context.tsv              every resistance/virulence hit in the controls,
                                     with coordinates and adjacent mobile elements
  pregate_per_genome.tsv           pre-gate call for each of the 25 controls
  pregate_hdc.tsv                  histidine decarboxylase screen, per genome
  truncation_controls.faa          synthetic controls used to test variant calling
```

`run_log.txt` records tool versions, retrieval dates, thresholds, the UniProt release the
reference set came from, the profile MD5 and self-score, the proximity window used for
the pre-gate, and every change made to the analysis after it began.

**The self-score matters.** Without it the percentage scale cannot be reconstructed, and
every number in the tables above becomes uncheckable. It is stored in three places:
`refs/self_score.txt`, `refs/gate2_profile.json`, and `run_log.txt`.

---

## Reproducing the analysis

```bash
bash   pipeline/00_setup.sh                     # tools + dbCAN databases, with checks
python pipeline/01_reference_proteins.py        # UniProt query → profile → self-score
python pipeline/02_build_panel.py --taxa taxa_list.txt --out panel/
bash   pipeline/03b_download_resume.sh data/screened_accessions.txt
bash   pipeline/pregate.sh                      # pre-gate: CARD/VFDB + HDC
python pipeline/pregate_context.py              # mobility context
python pipeline/04_validate_controls.py         # STOPS if controls fail
bash   pipeline/05_search_all.sh                # Gate 2
bash   pipeline/gate1_cazyme.sh                 # Gate 1
python pipeline/06_call_gates.py
python pipeline/07_summarize.py
python pipeline/test_truncation.py              # variant-calling controls
```

Step 04 is not optional. Without it a low observed prevalence cannot be distinguished
from a screen that does not work.

**`00_setup.sh` will stop if the dbCAN databases are incomplete.** Gate 1 needs all
three — the HMM profiles, the CAZy DIAMOND database, and dbCAN_sub. With only the HMM
database, `--tools all` silently degrades to a single method and the GH family counts
come out inflated. That is the error this revision corrects.

---

## Thresholds, fixed in advance

**Gate 2.** Each genome's highest full-sequence bit score against the profile is divided
by the profile self-score (1342.2 bits) and expressed as a percentage.

| Call | Criterion |
|---|---|
| **present** | normalized score > 50% |
| **absent** | normalized score ≤ 50% |
| *flagged* | normalized score 40–60% — called, but referred to enzyme assay |

Positive calls are sub-classified for reporting: **truncated** if the protein is shorter
than 80% of the model, **intact** at ≥ 90%, **point-variant** in between.

Ornithine, histidine and glutamate decarboxylases are **excluded from the reference set,
not included as decoys.** Whether a locus is present is not decided by which of several
similar profiles happens to score highest, and the margin between two related group II
PLP profiles carries no information about the target. Histidine decarboxylase is screened
separately at the pre-gate, where it belongs: in lactic acid bacteria it is
pyruvoyl-dependent and outside the group II family entirely.

**Gate 1.** A CAZyme family is counted only where at least two of dbCAN's three methods
agree. No pass threshold is applied — the gate returns a position on a distribution, not
a verdict. Any cut-off would be arbitrary, and qualification by the lactate/acetate route
is conditional on co-culture confirmation.

Species with zero detections are reported with a Wilson upper bound rather than as "0%".

---

## What the variant classification can and cannot do

`test_truncation.py` builds synthetic controls from a reference TyrDC — truncated to
95, 90, 85, 80, 70 and 50% of length, and with 5, 10, 20 and 40% of residues
substituted — and runs them through the calling rule.

| Control | Score | Call |
|---|---|---|
| intact | 95.9% | present / intact |
| 5% substituted | 90.1% | present / intact |
| 10–20% substituted | 73–83% | present / **point-variant** |
| 40% substituted | 46.2% | absent |
| truncated to 95% | 93.0% | present / intact |
| truncated to 80–90% | 79–88% | present / **point-variant** |
| truncated to 70% | 69.2% | present / **truncated** |
| truncated to 50% | 49.0% | absent |

All three sub-classifications fire, so a zero count in the real data is a real zero.

**The truncated class is only reached by severe truncation.** An allele cut to 80–90% of
model length scores high enough to be called, but the length ratio stays above 0.80 and
it is reported as a point-variant. Under the previous DIAMOND thresholds the truncated
class could not fire at all; it can now, but the boundary sits near 70% rather than at
the nominal 80%.

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

The pre-gate also screens for histidine decarboxylase, at a deliberately permissive
threshold (identity ≥ 40%, coverage ≥ 70%). A hit is a referral for culture confirmation,
not an exclusion — for a safety flag, missing a case is worse than a false alarm.

---

## Revision history

### August 2026 — method revision

Both gates were re-run after a co-author's review identified two problems.

**Gate 1 ran only one of dbCAN's three prediction methods.** `05_search_all.sh` passed
`--tools hmmer`, and `gate1_cazyme.sh` bypassed `run_dbcan` entirely and searched the HMM
database directly. Single-method domain hits overcount the repertoire: the *L. casei*
median was 26 under the old path and is 21 under the consensus rule. Every Gate 1 number,
Figure 2 and Table 2 were rebuilt.

**Gate 2 assigned each protein to whichever of several similar profiles scored highest.**
The reference set contained ornithine and histidine decarboxylases as decoys, neither of
which acts on levodopa. Presence is not decided by a ranking contest between related
profiles, and the margin between two group II PLP profiles says nothing about the target.
The reference set was rebuilt from UniProt (`reviewed:true`, `ec:4.1.1.*`,
`taxonomy_id:2`), a profile HMM built from the four sequences returned, and each genome
scored against the profile self-score.

Consequences: 570 → 633 detections, four → six carrier species, and the ambiguous
category disappeared. `APPLY.md` lists every file changed.

### During the original analysis

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

**The reference set is enterococcal-weighted.** Three of the four UniProt-reviewed
bacterial tyrosine decarboxylases are from *Enterococcus* and one from *L. brevis*; no
reviewed entry exists outside these two genera. Positive calls in more distant species
therefore score lower, and that distance is not by itself evidence of reduced activity.

*Lactiplantibacillus plantarum* IR BL0076 — the one isolate of this species reported to
carry *tyrDC* — has no public genome assembly and is absent from the panel, though the
locus itself is deposited as a sequence (GenBank JQ040309).

Detection of the locus is not evidence of activity. Functional confirmation requires
enzyme assay under controlled tyrosine availability.

---

## Citation

If you use this pipeline or these data, please cite the manuscript. Reference protein
accessions and the profile parameters are in Supplementary Table S5; the disclosure audit
of published Parkinson's disease trials is Supplementary Tables S3 and S4.

## Licence

Code: MIT. Data: CC-BY-4.0.
