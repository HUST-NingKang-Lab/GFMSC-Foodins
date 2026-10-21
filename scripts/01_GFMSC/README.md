# GFMSC construction scripts

These scripts support construction and annotation of the Global Food Microbial
Small Open Reading Frame Catalog (GFMSC).

## Scripts

1. `01_predict_short_ORFs.sh`: predict short proteins from MAGs or assembled
   metagenomic contigs with Pyrodigal.
2. `02_filter_smORFs.py`: retain ATG-initiated peptides of 5–50 aa after removing
   the initiator methionine and terminal stop symbol.
3. `03_build_GFMSC.sh`: combine the cFMD and fermented-grain peptide sets and
   cluster them with CD-HIT.
4. `04_assign_contig_taxonomy.sh`: assign taxonomy to assembled metagenomic
   contigs with MMseqs2 `easy-taxonomy` and GTDB.
5. `05_assess_smORF_novelty.sh`: search GFMSC representatives against sequence
   databases using the thresholds reported in the manuscript.

## Usage

Scripts accept input and output paths as command-line arguments and display
usage instructions when run without arguments. Run the scripts in numerical
order. Assembly inputs must be FASTA files; taxonomy and novelty searches require
prebuilt MMseqs2-compatible databases.

## Requirements

Python 3, Bash, Pyrodigal, CD-HIT and MMseqs2. Taxonomic assignment additionally
requires the GTDB database formatted for MMseqs2.
