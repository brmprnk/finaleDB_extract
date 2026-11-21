# finaleDB_extract: FinaleDB download and BAM reconstruction

This repository recreates, from public [FinaleDB](http://finaledb.research.cchmc.org/) records, the per-sample BAM files.
It downloads the hg38 fragment files of a list of FinaleDB records and rebuilds a sorted, indexed BAM from each with a modified
[Fragmentstein](https://github.com/uzh-dqbm-cmi/fragmentstein) (Balázs et al. 2024). No data are deposited here: the
records are public on FinaleDB, and the scripts and lists in this repository regenerate the BAMs deterministically.

Outputs of this repository are used for [comp](https://github.com/brmprnk/comp).

## Contents
1. [Scope](#1-scope)
2. [Requirements](#2-requirements)
3. [Setup](#3-setup)
4. [Walkthrough](#4-walkthrough)
5. [Input and output](#5-input-and-output)
6. [Runtime and disk](#6-runtime-and-disk)
7. [Records without hg38 data](#7-records-without-hg38-data)
8. [Provenance of the original BAMs](#8-provenance-of-the-original-bams)
9. [chr21 test set](#9-chr21-test-set)
10. [Repository layout](#10-repository-layout)
11. [Citation and license](#11-citation-and-license)

## 1. Scope

| Step | Here | Script |
|---|---|---|
| Record IDs of both cohorts | yes | `lists/cris_ids.txt`, `lists/jiang_ids.txt` |
| FinaleDB metadata (optional) | yes | `scripts/fetch_metadata.py` |
| Download the hg38 fragment files | yes | `scripts/download.sh` |
| Fragment file → sorted, indexed BAM (modified Fragmentstein) | yes | `scripts/convert.sh`, `fragmentstein/fragmentstein_finaledb.sh` |
| GC-bias correction, copy number, tumour fraction | no: downstream, see step 7 of the walkthrough | |

The ID lists cover two cohorts, Cristiano et al. 2019 (EGA EGAS00001003611; 534 FinaleDB records) and Jiang et al. 2015
(225 records). Any other list of FinaleDB records works the same way.

## 2. Requirements

All versions are pinned in `pixi.toml` and locked in `pixi.lock` (linux-64, osx-64, osx-arm64; bioconda has no Windows
builds of samtools or bedtools):

| Tool | Version | Used for |
|---|---|---|
| samtools | 1.21 | SAM/BAM conversion, sorting, indexing |
| htslib (bgzip, tabix) | 1.21 | the chr21 test set |
| bedtools | 2.31.1 | `bedToBam`, `bedtools getfasta` |
| awk | any POSIX awk; the environment adds GNU awk 5.3 | Fragmentstein's text processing (the scripts call `awk` from `PATH`) |
| curl | locked | downloads |
| Python | 3.12, with requests 2.32 | `scripts/fetch_metadata.py` |
| Fragmentstein | v2023.4a, modified copy in `fragmentstein/` | conversion; see [`fragmentstein/README.md`](fragmentstein/README.md) |
| deepTools (optional, `bigwig` environment) | 3.5 | bigWig coverage tracks |
| Reference genome | UCSC hg38 analysis set, 195 sequences | read sequences |

The original BAMs (Section 8) were made with the cluster's modules, without a container; their samtools and bedtools versions
are recorded in each BAM's `@PG` lines (Section 8). The Slurm wrappers run every step through `pixi run`, through the
Apptainer image built from `container.def` (the same locked environment), or with tools already on `PATH` (Section 3).

## 3. Setup

```bash
git clone https://github.com/brmprnk/finaleDB_extract && cd finaleDB_extract
pixi install                                   # the locked environment (https://pixi.sh)

# reference genome: UCSC hg38 analysis set (no alt contigs, decoys and EBV included, 195 sequences)
mkdir -p reference_genome
curl -o reference_genome/hg38.analysisSet.fa.gz \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
md5sum reference_genome/hg38.analysisSet.fa.gz  # 6d3c82e1e12b127d526395294526b9c8 (UCSC md5sum.txt)
gunzip reference_genome/hg38.analysisSet.fa.gz
pixi run samtools faidx reference_genome/hg38.analysisSet.fa
```

The original BAMs were made with this UCSC file: our copy's `hg38.analysisSet.2bit` has UCSC's md5 (`313b1f3476b471362748aceab1bf9991`).
Fragmentstein converts the sequence to upper case, so soft-masking does not matter.

On storage that does not allow pixi to set file permissions (SMB/CIFS mounts; `pixi install` fails with
"Operation not permitted"), use the container instead. It holds the same locked environment in one file:

```bash
apptainer build preprocessing.sif container.def
apptainer exec --bind /path/to/storage preprocessing.sif samtools faidx reference_genome/hg38.analysisSet.fa
```

For the Slurm wrappers, copy `slurm/site.env.example` to `slurm/site.env` (git-ignored) and fill in the scheduler
account, partition and QoS, `PREPROCESSING_DIR` (where the BAMs are written), `FRAGMENT_DIR`,
`GENOME`, and how the tools are run (`RUN="pixi run --frozen"`, `RUN="apptainer exec --bind /path/to/storage
preprocessing.sif"`, or empty with `module load` lines in the same file).
Load it in the submitting shell (`. slurm/site.env`) and submit from the repository root. The wrappers carry no
site-specific `#SBATCH` lines: sbatch takes the account, partition and QoS from `SBATCH_ACCOUNT`, `SBATCH_PARTITION`
and `SBATCH_QOS`.

## 4. Walkthrough

Shown for the Cristiano cohort; for Jiang replace `cris` by `jiang` and drop `-s` (step 4). Each step has a direct
command and a Slurm wrapper. Every step skips the records that are already complete, so it can be resubmitted after
an interruption.

**1. Record IDs.** `lists/cris_ids.txt` (534 records) and `lists/jiang_ids.txt` (225 records, one per sample) list
the FinaleDB records of the two cohorts. Every script also accepts any plain ID list (one `EE…` ID per line) or any
CSV with an `ID` column, such as a sample sheet.

**2. Metadata (optional).** Sample name, disease, read length and fragment count of each record from the FinaleDB API.
Records without an API result are flagged (Section 7):

```bash
pixi run python scripts/fetch_metadata.py lists/cris_ids.txt -o cris_metadata.csv
```

**3. Download** the fragment files (`<ID>.hg38.frag.tsv.bgz`, BED-like: chrom, start, end, MAPQ, strand):

```bash
pixi run bash scripts/download.sh -i lists/cris_ids.txt -o "$FRAGMENT_DIR/cris_frag" -j 4
sbatch slurm/download.sbatch cris                       # same, on the cluster
```

A file counts as complete when its size equals the server's; partial files are resumed. Records that FinaleDB serves
empty are listed in `empty_records.txt` and not kept (exit status 2).

**4. Convert** each fragment file to a sorted, indexed BAM with the modified Fragmentstein (read length 101, the
MAPQ filter of Section 8 per record, sequences from the reference):

```bash
pixi run bash scripts/convert.sh -i lists/cris_ids.txt -f "$FRAGMENT_DIR/cris_frag" -o "$PREPROCESSING_DIR/cris" \
    -g reference_genome/hg38.analysisSet.fa -s lists/cris_conversion_settings.csv -j 16 -t 4
sbatch slurm/convert.sbatch cris                        # same: 64 CPUs, 16 records at a time
```

`-s lists/cris_conversion_settings.csv` reproduces the per-record MAPQ filter of the original Cristiano BAMs
(Section 8). Without `-s` every record keeps fragments with MAPQ > 5, as for all Jiang records. One log per record is
written to `logs/convert/` (`logs/convert_<cohort>/` for the Slurm wrapper), and failed records are listed in
`failed_records.txt` there.

**5. Index.** Step 4 indexes every BAM. `scripts/index.sh -d DIR` (`sbatch slurm/index.sbatch cris`) only repairs a
folder, for example BAMs copied without their `.bai`.

**6. Outputs.** `<PREPROCESSING_DIR>/cris/<ID>.hg38.frag.tsv.bam` and `.bam.bai`, and the same under
`<PREPROCESSING_DIR>/jiang/` (Section 5). To check that a rebuilt BAM has the same reads as another copy (headers and
compression ignored): `pixi run bash scripts/compare_bams.sh A.bam B.bam`.

**7. Next steps** (not included here). The BAMs work with the usual cfDNA tools, for example GC-bias correction
with GCfix, or copy number and tumour fraction with HMMcopy `readCounter` and ichorCNA.

## 5. Input and output

- **Input:** FinaleDB record IDs (`EE` + number), as a plain list or the `ID` column of a CSV.
- **Output:** `<PREPROCESSING_DIR>/<cohort>/<ID>.hg38.frag.tsv.bam` (with the Slurm wrappers `<cohort>` is `cris` or
  `jiang`), each with its index `<ID>.hg38.frag.tsv.bam.bai`.
- **BAM contents:** coordinate-sorted; 195 `@SQ` lines of hg38.analysisSet; two reads per fragment (flags 99 and 147),
  each `min(101, fragment length)` bases; MAPQ = the FinaleDB MAPQ; TLEN = the fragment length; base qualities `F`.

## 6. Runtime and disk

| Per record (about 40 million fragments) | Size |
|---|---|
| Fragment file (download) | 0.20–0.28 GB |
| BAM + index | median 2.45 GB (Cristiano, range 0.97–8.58), 2.13 GB (Jiang, range 1.37–4.17) |
| Temporary BED and SAM files during the conversion | roughly 12 GB (estimate, not measured) |
| All BAMs | 1.44 TB (Cristiano, 534 records), 0.49 TB (Jiang, 225 records) |

Runtime: in the 2025 runs on one 64-CPU node (64 records in parallel, 4 samtools threads each, bigWig included), the
225 Jiang BAMs were all finished within 5.2 h of the first one, and 195 Cristiano BAMs within 5.7 h. That is about
35–45 records per node-hour, or 1–1.5 h per record under that load. The download is limited by the FinaleDB server
(about 0.25 GB per record). `scripts/convert.sh` logs the minutes each record takes.

## 7. Records without hg38 data

The Jiang et al. 2015 study has 249 FinaleDB records, EE85754–EE86002. For 24 of them FinaleDB has no hg38 data: the
API returns an empty result and the hg38 fragment file is a 28-byte BGZF end-of-file block with no fragments. The 24
records are EE85778, EE85782, EE85801, EE85803, EE85804, EE85819, EE85820, EE85824, EE85850, EE85852, EE85889,
EE85894, EE85909, EE85910, EE85914, EE85929, EE85934, EE85954, EE85958, EE85976, EE85977, EE85979, EE85993 and
EE85995. Downloads in March 2026 and checks in November 2026 returned the same empty files, so this is not a transfer error.

The 2025 run converted all 249 records, and the empty files gave 24 BAMs without reads (mean coverage 0). Those BAMs
have since been deleted, and the 24 records are not in `lists/jiang_ids.txt`. The metadata query of that run also had a bug: for an empty API result it repeated
the previous record's row, so the 24 records appeared as 24 duplicate rows of their neighbours: 249 rows for 225 records, one record per
sample. The 225 records of `lists/jiang_ids.txt` all have data.

Both checks are now part of the normal path. `scripts/download.sh` detects empty records, lists them in
`empty_records.txt` and keeps no file. `fragmentstein_finaledb.sh` refuses an empty input (exit status 3).
`scripts/fetch_metadata.py` writes one row per requested ID, with `status=no_result` for these records. No separate
re-download step is needed.

## 8. Provenance of the original BAMs

The original BAMs are the 759 BAMs of the two ID lists that this repository was written to reproduce; they were
used for comp. Each BAM's `@PG` line records how it was made:

| Records | Fragmentstein copy | MAPQ filter | samtools, bedtools | When |
|---|---|---|---|---|
| Cristiano, 262 records | `stein.sh` | 0 (keeps MAPQ ≥ 1) | 1.16.1–1.20, 2.31.1 (one record 2.30.0) | Aug–Oct 2024 |
| Cristiano, 272 records | `stein_parallel.sh` | 5 (keeps MAPQ ≥ 6) | 1.21, 2.31.1 | Aug–Oct 2025 |
| Jiang, all 225 records | `stein_parallel.sh` | 5 (keeps MAPQ ≥ 6) | 1.21, 2.31.1 | Apr 2025 |

All used read length 101, sorting on and the same reference. The exact differences of both copies from upstream
Fragmentstein are in `fragmentstein/history/`, and `fragmentstein/README.md` describes them.
`lists/cris_conversion_settings.csv` holds the MAPQ filter of each Cristiano record, as read from its BAM, so step 4
rebuilds the original reads. samtools versions change the compression and the header, not the reads.

## 9. chr21 test set

A small end-to-end test, and a small data set with both classes for testing downstream analyses.

- **Records:** `lists/chr21_test.csv`, 20 Cristiano records: 10 healthy, 5 breast cancer and 5 colorectal cancer.
  In each class the records are
  those whose fragment count is closest to the class median: 37.6–38.1 million for the healthy records, 41.8–42.9
  million for breast cancer and 37.2–40.5 million for colorectal cancer.
- **What it does:** fetches only the chr21 fragments of each record from FinaleDB with `tabix`, which transfers
  about 3.5 MB per record instead of the 0.25 GB file. It converts them with the same Fragmentstein copy and MAPQ
  settings, indexes the BAMs, and writes `sizes.tsv` and `checksums.md5` (the md5 of each BAM's reads, header
  excluded, to compare two builds).
- **Size:** measured from the indexes of the original BAMs, chr21 holds 1.27 % of the reads. That is 3.5–3.9 MB of
  fragments (about 550,000) and 29–32 MB of BAM per record, about 75 MB of fragment files and 600 MB of BAMs for all
  20 records. These are not hosted anywhere. The script rebuilds them from FinaleDB: on the cluster each record took under a
  minute (4 records in parallel, 2 threads each), and the build came to 66.6 MB of fragment files and 599.8 MB of
  BAMs.

```bash
pixi run bash scripts/make_chr21_test.sh -g reference_genome/hg38.analysisSet.fa -o chr21_test -j 4
sbatch slurm/chr21_test.sbatch                          # same, output in <PREPROCESSING_DIR>/chr21_test
```

The output folder has the same layout as `PREPROCESSING_DIR` (`chr21_test/cris/<ID>.hg38.frag.tsv.bam`), so it
can stand in for it in downstream tests. Read names number the fragments of the chr21 file, so they differ from those of the same
reads in a full BAM. All other fields are the same.

## 10. Repository layout

```
scripts/                 download.sh, convert.sh, index.sh, compare_bams.sh, make_chr21_test.sh, fetch_metadata.py
fragmentstein/           the modified Fragmentstein, its diff from upstream and the history of the original copies
lists/                   record IDs (cris, jiang, chr21 test) and the per-record MAPQ filter of the original Cristiano BAMs
slurm/                   sbatch wrappers; site settings in slurm/site.env (template slurm/site.env.example)
pixi.toml, pixi.lock     the locked environment
container.def            Apptainer image with that environment
```

## 11. Citation and license

When you use this repository, please cite it (`CITATION.cff`) and Fragmentstein: Balázs Z, Gitchev T, Ivanković I,
Krauthammer M. Fragmentstein — facilitating data reuse for cell-free DNA fragment analysis. *Bioinformatics*
40(1):btae017 (2024), doi:[10.1093/bioinformatics/btae017](https://doi.org/10.1093/bioinformatics/btae017).
The data are those of FinaleDB (Zheng et al., *Bioinformatics* 37(16):2502–2503, 2021, doi:[10.1093/bioinformatics/btaa999](https://doi.org/10.1093/bioinformatics/btaa999))
and of the original studies: Cristiano et al., *Nature* 2019, doi:[10.1038/s41586-019-1272-6](https://doi.org/10.1038/s41586-019-1272-6),
and Jiang et al., *PNAS* 2015, doi:[10.1073/pnas.1500076112](https://doi.org/10.1073/pnas.1500076112).

Distributed under the GNU General Public License v3.0 ([LICENSE](LICENSE)), as is upstream Fragmentstein.
