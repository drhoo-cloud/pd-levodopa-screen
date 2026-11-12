# pd-levodopa-screen

Data and code behind *Genome screening for levodopa compatibility and
gut–brain axis support before probiotic selection in Parkinson's disease*,
submitted to *npj Biofilms and Microbiomes*.

Genome sequences are **not** deposited here; every assembly is public in
RefSeq and recoverable from the accession columns of the per-genome files.

## Verify everything with one command

```
python3 verify_all.py
```

It recomputes **53 of the numbers printed in the paper** directly from the
files in `data/` and prints, for each, the value stated in the paper beside
the value it just calculated. Python 3 only — no third-party packages, no
network. It exits non-zero if anything disagrees.

## What changed in v3.0.0

The screen was extended from 200 genomes per species to **every genome
deposited in RefSeq for the 22 species**: 19,771 at Gate 2 and 19,765 at
Gate 1, against 5,372 before. Files that the census superseded
(`per_genome_calls.tsv`, `gate1_gh_families.txt`, `table_by_species.tsv`,
`screened_accessions.txt`) are gone from this release and remain available
at the v2.1.0 DOI.

`refs/` and `pipeline/` are unchanged from v2.1.0. The profile that produced
these calls is still `refs/tyrdc.hmm`, md5 `e5ddd734f368…`, self-score 1342.2.

## What is where

| Claim in the paper | File |
|---|---|
| 19,771 genomes at Gate 2; six carrier species; five with the operon | `data/gate2_calls_per_genome.tsv`, `data/gate2_by_species.tsv` |
| 19,765 genomes at Gate 1; medians 7–28 glycoside hydrolase families | `data/gate1_gh_families_per_genome.tsv` |
| Operon context of every carrier | `data/locus_context_census.tsv`, `data/locus_context_by_species.tsv` |
| Polysaccharide lyases absent from twenty of twenty-two species | `data/cazyme_all_categories.tsv` |
| *L. plantarum* 0 / 1,665, decarboxylase family present in 98.62% | `data/lp_decarboxylase_census.tsv`, `data/lp_pf00282_in_prophage.tsv` |
| Prophages in 98.3% of genomes, 7,305 regions, 27 family members inside one | `data/lp_prophage_census.tsv`, `data/lp_pf00282_in_prophage.tsv` |
| 1,665 genomes collapse to 161 lineages; bound widens to 2.33% | `data/lp_dereplication_census.tsv`, `data/lp_cluster_membership.tsv`, `scripts/cluster_single_linkage.py` |
| Prevalence by isolation source, 1.28% to 34.53% | `data/isolation_source.tsv`, `data/sample_metadata.tsv` |
| Profile detects 145 of 147; four of 211 outside Bacillota; 50.00% excluded | `data/profile_validation_scores.tsv` |
| 166 related decarboxylases, none called, highest 4.72% | `data/profile_benchmark_scores.tsv` |
| Pre-gate on 25 control genomes, 123 *van* hits, 45.8% background | `data/pregate_context.tsv`, `data/pregate_per_genome.tsv`, `scripts/verify_pregate.py` |
| Neither gate predicts the other | `data/twogate_census_stats.tsv`, `scripts/recompute_twogate_census.py` |
| The 22 trial components resolved to strain | `data/disclosure_audit_named_strains.tsv` |

## Reproducing a Gate 2 call yourself

```
hmmsearch -E 1e-5 refs/tyrdc.hmm <your_proteome.faa>
```

A hit counts as present when its bit score exceeds **half the profile's
self-score** (`refs/self_score.txt`, 1342.2), i.e. above 671.1. The canonical
operon additionally requires *tyrS* and *tyrP* within 20 kb.
`refs/gate2_profile.json` states the model and every threshold.

## Notes a reader should have

- `data/gate2_calls_per_genome.tsv` holds **19,772 rows**. One of them,
  *Enterococcus faecalis* V583 (`GCF_000007785.1`), is the control genome the
  self-score is defined against and is **excluded from the census**; the
  `in_census` column marks it. Filtering `in_census == yes` gives the 19,771
  of the paper.
- Gate 1 covers 19,765 genomes and Gate 2 covers 19,771. The two gates were
  run on separate retrievals and six assemblies reached one gate but not the
  other. Both figures are stated in the paper; neither is a subset error.
- Dereplication at 99.5% ANI was computed on the full 1,665-genome census for
  *L. plantarum*, and on the capped set for the other species.
- `MANIFEST.tsv` lists every file with its line count, byte size and MD5.
- All text files use LF line endings and UTF-8.

## Layout

```
data/       per-genome and per-species results
refs/       the profile, its reference sequences, thresholds and the two trees
pipeline/   the pipeline as it was run for the original screen
scripts/    the analyses that are not a single command
verify_all.py   recomputes the paper's numbers from the files above
run_log.txt     the order the analysis was run in, with the date of each step
run_log_laboratory_record.txt
                the laboratory's working log, translated from Korean
MANIFEST.tsv    line count, size and MD5 of every file
CHANGELOG.md    what changed in this release
```

## Licence

MIT, see `LICENSE`.
