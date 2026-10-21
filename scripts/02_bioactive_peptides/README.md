# Bioactive peptide prediction scripts

These scripts support Foodin prediction, prediction of other bioactive peptide
classes and comparisons with public peptide resources.

## Scripts

1. `01_run_APEX.sh`: run APEX on the non-redundant GFMSC representatives.
2. `02_run_ProTrek_similarity.py`: encode GFMSC and reference peptides with
   ProTrek-650M and retain sequence pairs above the selected similarity score.
3. `03_extract_ProTrek_candidates.py`: export unique GFMSC candidates from the
   retained ProTrek sequence pairs.
4. `04_merge_activity_fastas.py`: merge AMP, ACP, AOP, AIP, ADP and AHP FASTA
   files into one activity table.
5. `05_compare_candidates_to_references.sh`: compare candidates with public
   reference collections using MMseqs2.
6. `06_map_foodins_to_long_ORFs.sh`: map Foodins to ORFs of at least 52 aa and
   export alignment coordinates and coverage.
7. `07_select_best_hit.awk`: select the lowest-E-value match for each query.

## Usage

Scripts accept input and output paths as command-line arguments. Run APEX once
for the GFMSC representatives and run the ProTrek workflow separately for each
non-antimicrobial activity class before merging the final FASTA files.

## Requirements

Python 3, Bash, AWK, APEX, ProTrek-650M, PyTorch, tqdm and MMseqs2. APEX and
ProTrek model weights must be obtained from their respective projects.
