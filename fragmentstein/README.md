# Modified Fragmentstein

`fragmentstein_finaledb.sh` converts one FinaleDB hg38 fragment file to a sorted, indexed BAM. It is a modified copy
of **Fragmentstein v2023.4a** (`scripts/fragmentstein.sh` at upstream commit
[`d62329c`](https://github.com/uzh-dqbm-cmi/fragmentstein/commit/d62329cffc6e38d92b3bf7ad0a9ae6012a82eb3b),
18 July 2023, the last change to that script) by Zsolt Balázs and Todor Gitchev,
<https://github.com/uzh-dqbm-cmi/fragmentstein>.

- **License:** upstream is GPL-3.0. This copy is a derivative work and is distributed under GPL-3.0 too, like the
  rest of this repository (`LICENSE`).
- **Citation:** Balázs Z, Gitchev T, Ivanković I, Krauthammer M. Fragmentstein — facilitating data reuse for cell-free
  DNA fragment analysis. *Bioinformatics* 40(1):btae017 (2024). doi:[10.1093/bioinformatics/btae017](https://doi.org/10.1093/bioinformatics/btae017)

## What Fragmentstein does

Every FinaleDB fragment (chrom, start, end, MAPQ, strand) with MAPQ above the filter becomes a properly paired
read pair: a forward read (flag 99) at the fragment start and a reverse read (flag 147) ending at the fragment end,
each `min(101, fragment length)` bases long. The sequences come from the reference genome, the base qualities are
all `F` (Phred 37), the MAPQ is the FinaleDB MAPQ, and the template length is the fragment length. Read names are
`F1`, `F2`, ... in the order of the fragment file.

## Changes with respect to upstream

The complete diff is [`upstream.diff`](upstream.diff). The conversion commands (the `awk`, `bedToBam`,
`bedtools getfasta` and `samtools` pipelines) are upstream's; the changes are around them:

| # | Change | Effect on the reads |
|---|---|---|
| 1 | Defaults set to the values of the original BAMs: read length 101, MAPQ filter 5 (`-qf`, keeps MAPQ > 5), default MAPQ 5 (`-qd`), base quality F, sort and index on. Upstream defaults: `-qf 30 -qd 60`, no sorting | the values used; `-qf` and `-qd` stay options |
| 2 | FinaleDB `.tsv.bgz` input only (BED/BEDPE branches removed); `-o` required | none |
| 3 | Temporary files `<output>.tmp.bed/.tmp.sam` without the `.bam` in the middle, optionally in `--temp DIR` | none |
| 4 | Optional deepTools `bamCoverage` bigWig (`--bigwig`) | none (extra file) |
| 5 | `@PG` line: `VN:2023.4a-finaledb-1.0.0` and the script name | header only |
| 6 | Checks: empty input stops (exit 3); the BAM must hold one read per line of the temporary BED; partial BAMs are removed on failure | none |

`-qd` does not change the reads: `bedToBam -mapq` sets a placeholder MAPQ that the next step replaces with the
FinaleDB MAPQ column. It does appear in the `@PG` line.

## Provenance of the published BAMs

The original BAMs (main README, Section 8) were made with two earlier batch copies of the same upstream script. Their exact differences from
upstream are in [`history/`](history/):

| Copy | Used for | Parameters (from each BAM's `@PG` line) | Tools |
|---|---|---|---|
| `stein.sh` ([diff](history/stein_2024.diff)), Aug–Oct 2024 | 262 Cristiano records (and one record outside `lists/cris_ids.txt`) | `-r 101 -qf 0 -qd 0 -s true -t 4` | samtools 1.16.1–1.20, bedtools 2.31.1 (one record 2.30.0) |
| `stein_parallel.sh` ([diff](history/stein_parallel_2025.diff)), Apr–Oct 2025 | the other 272 Cristiano records and all 225 Jiang records | `-r 101 -qf 5 -qd 5 -s true -t 4` | samtools 1.21, bedtools 2.31.1 |

Both copies hard-code the parameters and loop over a folder, `stein_parallel.sh` in parallel with `xargs`. Both also
write a deepTools bigWig, and `stein_parallel.sh` leaves `VN` and the script name in the `@PG` line empty, because
`VERSION` is not exported to the `xargs` subshells. Otherwise their conversion commands are upstream's.
`lists/cris_conversion_settings.csv` records the MAPQ filter of every Cristiano record, and `scripts/convert.sh -s`
applies it. With that list, `fragmentstein_finaledb.sh` should rebuild the published reads; check it with
`scripts/compare_bams.sh`. Different samtools versions change the BGZF compression and the `@PG` lines, not the reads.
