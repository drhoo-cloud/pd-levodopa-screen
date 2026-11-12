# v3.0.0 — the screen extended to every deposited genome

The screen was run on 200 genomes per species. This release runs it on
**every genome deposited in RefSeq for the 22 species**: 19,771 at Gate 2 and
19,765 at Gate 1, against 5,372 before.

No conclusion reversed. *Lactiplantibacillus plantarum* went from 0/1,521 to
**0/1,665**, which narrows its zero-prevalence upper bound from 0.25% to
0.23%. The pooled prevalence across the panel rose from 11.8% to 57.5%,
because the cap had been truncating the two *Enterococcus* species hardest —
which is the point of reporting the census rather than a sample.

## What is new

- **Gate 2 across the census.** `data/gate2_calls_per_genome.tsv` replaces
  `per_genome_calls.tsv`: 19,771 genomes plus the V583 control, marked by an
  `in_census` column. Six carrier species, five with the canonical operon.
- **Gate 1 across the census.** `data/gate1_gh_families_per_genome.tsv`
  replaces `gate1_gh_families.txt`: 19,765 genomes, species medians 7 to 28.
- **Operon context for every carrier**, not a subsample:
  `data/locus_context_census.tsv`, 11,365 rows.
- **Dereplication recomputed on the full census.** The earlier figure was
  computed on the 1,521 *L. plantarum* genomes available at the time of the
  first deposit, before 144 were added. All 1,665 were re-clustered with Mash
  at 99.5% ANI: **161 independent lineages**, upper bound 2.33%, replacing
  156 and 2.40%.
- **Table S4b recomputed on the census.** The association models had been run
  on 5,372 genomes and are now run on 19,765:
  `data/twogate_census_stats.tsv`, with `scripts/recompute_twogate_census.py`
  and both marker trees in `refs/`.
- **`verify_all.py`.** One command recomputes 53 of the numbers printed in
  the paper from the files in `data/` and reports any that disagree. Python 3
  only, no third-party packages, no network.
- **`MANIFEST.tsv`** lists every file with its line count, size and MD5.
- **`run_log.txt`** is now an English record of all eleven steps with dates.
  The laboratory's original working log, previously the Korean `run_log.txt`,
  is translated in `run_log_laboratory_record.txt`; the Korean original
  remains in this repository's history at v2.1.0.

## What was removed

`per_genome_calls.tsv`, `gate1_gh_families.txt`, `table_by_species.tsv` and
`screened_accessions.txt` are superseded by the census files and are not in
this release. They remain available at the v2.1.0 DOI.

## What is unchanged

`refs/` and `pipeline/` are byte-identical to v2.1.0. The profile that
produced these calls is still `refs/tyrdc.hmm`, md5 `e5ddd734f368…`,
self-score 1342.2 bits. Two files are added to `refs/`: the rpoB and 16S
species trees used for the phylogenetic models.

`data/pregate_context.tsv` and `data/pregate_per_genome.tsv` are unchanged
apart from line endings, which are now LF throughout the repository.
