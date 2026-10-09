# Bioactive peptide prediction scripts

This directory contains prediction and sequence-comparison steps that
support the current manuscript.

## Scripts

1. `01_run_APEX.sh`: run APEX on the non-redundant GFMSC representatives.
2. `02_run_ProTrek_similarity.py`: encode GFMSC and reference peptides with
   ProTrek-650M, calculate pairwise similarity scores and retain matches at the
   selected score threshold. Run this script separately for ACP, AOP, AIP, ADP
   and AHP reference FASTA files.
3. `03_extract_ProTrek_candidates.py`: reduce the retained sequence pairs to
   unique GFMSC candidates and export the corresponding activity-specific
   FASTA file.
4. `04_merge_activity_fastas.py`: merge the final AMP, ACP, AOP, AIP, ADP and
   AHP FASTA files into one table, retaining multiple predicted activities in
   one row per peptide.
5. `05_compare_candidates_to_references.sh`: compare a candidate peptide set
   with a corresponding public reference collection using MMseqs2.
6. `06_map_foodins_to_long_ORFs.sh`: map Foodins to the collection of ORFs of
   at least 52 aa and retain alignment coordinates and coverage values.
7. `07_select_best_hit.awk`: select the lowest-E-value match for each query
   from the extended MMseqs2 output.

