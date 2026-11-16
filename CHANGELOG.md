# v3.2.0 — the isolation-source question, answered as a negative

Referee 2 asked for prevalence stratified by isolation source. Pooled across
the census it looks like a large effect, and it is not one: it is species
composition. This release deposits the per-genome source calls, the
within-species tests and the model, so the point can be checked either way.

- **`data/census_isolation_source.tsv`** — isolation source for all
  **19,771** census genomes, taken from the BioSample record
  (`isolation_source`, falling back to `host`) and classified under the same
  rules as the 5,368-genome `isolation_source.tsv`.
- **`data/census_source_by_species.tsv`**, **`census_source_by_species_long.tsv`**
  — the within-species tests reported in Supplementary Table S4c, and the
  full species x source cross-tabulation behind them.
- **`data/copynum_by_species.tsv`** — tyrDC copy number per species. This
  backs the multi-copy column of Supplementary Table S3a, which had no
  deposited source until now.
- **`scripts/source_model.py`**, **`scripts/source_by_species.py`** — the
  logistic models and the exact per-species tests.

## What the census shows

| Isolation source | Genomes | Carriers | Prevalence | *Enterococcus* share |
|---|---|---|---|---|
| probiotic / product | 132 | 16 | 12.12% | 10.6% |
| food / fermented | 2,356 | 533 | 22.62% | 19.3% |
| unknown / not stated | 4,644 | 2,539 | 54.67% | 52.4% |
| human / clinical | 7,736 | 4,755 | 61.47% | 61.2% |
| other (stated) | 3,435 | 2,414 | 70.28% | 69.1% |
| animal | 835 | 625 | 74.85% | 72.6% |
| environment / industrial | 633 | 483 | 76.30% | 72.0% |

Prevalence tracks the *Enterococcus* share of the category rather than the
source: across the seven categories the two agree at **r = 0.999**. Under the
specification the paper uses for the two gates in Supplementary Table S4b — a
logistic model with species fixed effects, fitted to the species containing
both outcomes — source adds nothing: **likelihood ratio chi2 = 7.69, df = 6,
P = 0.26, n = 12,277**, and no individual source term approaches
significance (lowest P = 0.085). Without the species term the same comparison
on the same data gives **chi2 = 1,820.3**, which measures the confounding and
not an effect.

## Within a species

Several species x source cells hold fewer than five genomes, where the
chi-square approximation is not safe. The per-species P value is therefore
exact: the source labels are permuted within the species, 50,000 times, seed
20261112, and Holm-corrected across the six species that contain both
outcomes. The range and spread are restricted to the sources holding at
least ten genomes of that species, because a category holding one genome
returns 0% or 100% and that is not a spread.

| Species | Genomes | Range across sources | Spread | Exact P | Holm P |
|---|---|---|---|---|---|
| *Latilactobacillus curvatus* | 118 | 12.50% – 71.43% | 58.93 pp | 0.003 | 0.019 |
| *Levilactobacillus brevis* | 306 | 68.24% – 96.43% | 28.19 pp | 0.011 | 0.053 |
| *Enterococcus faecalis* | 5,076 | 97.44% – 100.00% | 2.56 pp | 0.040 | 0.162 |
| *Limosilactobacillus reuteri* | 533 | 7.35% – 17.14% | 9.79 pp | 0.271 | 0.814 |
| *Enterococcus faecium* | 6,000 | 99.64% – 100.00% | 0.36 pp | 0.800 | 1.000 |
| *Pediococcus pentosaceus* | 244 | 0.00% – 2.91% | 2.91 pp | 0.889 | 1.000 |

One of the six survives correction. The asymptotic test would have called
three of them significant, *Enterococcus faecalis* among them, where the
entire spread is 2.56 percentage points on 5,076 genomes. The species that
does survive, *Latilactobacillus curvatus*, rests on 3 of 24 human isolates
against 10 of 14 in the residual category, from 118 genomes in all.

Pooling across species is what the paper declines elsewhere for the
panel-wide prevalence, so the main text reports no by-source prevalence;
Supplementary Table S4c carries the within-species figures instead.

`verify_all.py` checks these too: **80 numbers reproduced, 0 differing.**

Nothing from v3.1.0 is changed or removed. One thing is retired: the
per-species rows that `data/source_model_stats.tsv` carried in preparation
for this release used the asymptotic test, which these cell counts do not
support. `census_source_by_species.tsv` replaces them.

# v3.1.0 — the alternative route deposited, and two wording corrections

The tyrDC gate sees one route. A distinct levodopa-metabolising pathway runs
through commensal *Bifidobacterium*, by aromatic aminotransferase and then
aryl-lactate dehydrogenase, and the profile does not detect it. That census
was run for the referee response but never deposited. It is here now.

- **`data/bifido_aldh_aat.tsv`** — per-genome aminotransferase and
  dehydrogenase calls for all **1,872** genomes of the four panel
  *Bifidobacterium* species, with the distance between the two genes and a
  20 kb adjacency call.
- **`data/bifido_by_species.tsv`** — the same census per species.

The aminotransferase is near-universal (**1,868 of 1,872**), so the
discriminating gene is the dehydrogenase, present in **1,623 of 1,872**
(86.7%): 920 of 922 *B. longum*, 339 of 340 *B. breve*, 364 of 365
*B. bifidum*, and **none of 245 *B. animalis***. That split reproduces the
reported culture phenotype in all four species.

`verify_all.py` now checks these too: **59 numbers reproduced, 0 differing.**

## Two wording corrections

Neither changes a number.

- The v3.0.0 note said the screen "was run on 200 genomes per species". The
  200 cap applied to 17 of the 22 species. Four had fewer than 200 genomes in
  RefSeq, and *Lactiplantibacillus plantarum* was run at full depth from the
  start — which is why the v3.0.0 note could say in the same breath that it
  went from 0/1,521. `CHANGELOG.md` and `README.md` now say "at most 200 per
  species, except *L. plantarum*".
- Two file names were wrong in the manuscript's referee response, not here,
  and are noted for anyone following a citation: the disclosure audit is
  `data/disclosure_audit_named_strains.tsv` (not `disclosure_audit_all27.tsv`)
  and the pre-gate context table is `data/pregate_context.tsv` (not
  `pre_gate_context_table.tsv`).

Nothing in `data/` from v3.0.0 is changed or removed, and `refs/` and
`pipeline/` remain byte-identical to v2.1.0.

# v3.0.0 — the screen extended to every deposited genome

The screen was run on at most 200 genomes per species, except
*Lactiplantibacillus plantarum*, which was run at full depth from the start.
This release runs it on **every genome deposited in RefSeq for the 22
species**: 19,771 at Gate 2 and 19,765 at Gate 1, against 5,372 before.

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
