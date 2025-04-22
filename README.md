# finaleDB_extract
Download anonymized cfDNA files in frag.tsv.bgz file format and important sample metadata.
Then convert those to .bam files using Fragmentstein in a parallelized manner, including coverage .bw files to assert coverage tracks.

## Requirements
- Python packages ```requests, pandas, bioconda::deeptools```
- A reference genome (hg38). The code by default expectes a folder /reference_genome in the project directory containing hg38.analysisSet.fa, hg38.analysisSet.fa.fai as reference genome files.
These file names can be edited in stein.sh or stein_parallel.sh for that matter.

## Instructions
1.
Run finale_db.py with the range of sample ID's specified and the name of the output csv specified.
This retrievesa a sample_info.csv where all the sample information is summarized.

2.
Run download_ee.sh, where the sample_info.csv is given to download all the included sample_ids.
This can take ~250Mb per file, so beware of long downloading times and large storage requirements.

3.
Run stein.sh or stein_parallel.sh to map the donwloaded tsv files to bam files with the same name.
They contain the same logic but stein_parallel supports multiple processes (note CPU requirements and larger memory requirements).

# Author
Bram Pronk -- I.B.Pronk@tudelft.nl

## Citations

[1] Haizi Zheng, Michelle S Zhu, Yaping Liu, FinaleDB: a browser and database of cell-free DNA fragmentation patterns, Bioinformatics, Volume 37, Issue 16, August 2021, Pages 2502–2503, https://doi.org/10.1093/bioinformatics/btaa999

[2] Zsolt Balázs, Todor Gitchev, Ivna Ivanković, Michael Krauthammer, Fragmentstein—facilitating data reuse for cell-free DNA fragment analysis, Bioinformatics, Volume 40, Issue 1, January 2024, btae017, https://doi.org/10.1093/bioinformatics/btae017
