# GFMSC construction scripts

These scripts support the construction and annotation of the Global Food Microbial Small Open Reading Frame Catalog (GFMSC).

## Scripts

1. `01_predict_short_ORFs.sh`: predict short proteins from MAGs or assembled
   metagenomic contigs with Pyrodigal.
2. `02_filter_smORFs.py`: retain ATG-initiated peptides of 5-50 aa after
   removal of the initiator methionine and terminal stop symbol.
3. `03_build_GFMSC.sh`: combine the cFMD and fermented-grain peptide sets and
   cluster them with CD-HIT.
4. `04_assign_contig_taxonomy.sh`: assign taxonomy to assembled metagenomic
   contigs with the MMseqs2 `easy-taxonomy` workflow and GTDB.
5. `05_assess_smORF_novelty.sh`: search GFMSC representatives against one or
   more sequence databases using the thresholds reported in the manuscript.
