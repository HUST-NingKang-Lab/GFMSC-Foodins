# Global Food Microbial SmORF Catalog (GFMSC)

This directory contains scripts for short open reading frame (smORF) prediction,
peptide filtering, catalog construction, contig-level taxonomic assignment and
sequence homology searches for the Global Food Microbial SmORF Catalog (GFMSC).
Bioactivity prediction and downstream peptide analyses are organized separately
in `02_bioactive_peptides/`.

## Input data

The study combines 10,112 prokaryotic metagenome-assembled genomes (MAGs) from
curatedFoodMetagenomicData (cFMD) and 322 fermented-grain metagenomes deposited
in the Genome Sequence Archive under accession **CRA028755**.

This workflow starts from curated MAG assemblies and host-depleted metagenomic
assemblies. Read quality control, host removal, assembly and MAG generation are
upstream steps and are not performed by the scripts in this directory. Complete
assemblies, reference databases and source MAG metadata must be obtained
separately.

## Scripts

| Script | Purpose | Main output |
| --- | --- | --- |
| `01_predict_short_ORFs.sh` | Predict proteins from nucleotide assemblies with short-gene settings | One `<source>.proteins.faa` file per input assembly |
| `02_filter_smORFs.py` | Process predicted proteins and retain peptides of 5–50 amino acids | Merged, filtered peptide FASTA for each dataset |
| `03_build_GFMSC.sh` | Combine cFMD and fermented-grain peptides and cluster them with CD-HIT | Redundant peptides, family representatives and `.clstr` membership file |
| `04_assign_contig_taxonomy.sh` | Assign taxonomy to metagenomic contigs using MMseqs2 and GTDB | Contig-level taxonomic assignments and reports |
| `05_assess_smORF_novelty.sh` | Search representative peptides against reference sequence databases | One alignment table per database |

## Requirements

Run on Linux with Bash, GNU `find` and GNU `xargs`. All software executables must
be available in `PATH`. The filtering script uses the Python standard library;
no additional Python packages are required for that script.

| Software | Version reference | Executable |
| --- | --- | --- |
| Python | 3.9 or newer recommended | `python3` |
| Pyrodigal | 3.6.3 recorded in the analysis notes | `pyrodigal` |
| CD-HIT | 4.8.1 | `cd-hit` |
| MMseqs2 | 18.8cc5c | `mmseqs` |

## Running the workflow

The examples below assume this directory is located at the repository root as
`01_GFMSC/`. Run commands from the repository root. Input and output paths are
illustrative and can be replaced with local paths. Use fresh output directories
for independent analyses.

### 1. Predict short ORFs

```bash
bash 01_GFMSC/01_predict_short_ORFs.sh \
  data/cFMD_MAGs results/01_GFMSC/proteins/cFMD 8

bash 01_GFMSC/01_predict_short_ORFs.sh \
  data/fermented_grain_assemblies results/01_GFMSC/proteins/grain 8
```

The script recursively finds `.fa`, `.fna` and `.fasta` files and calls Pyrodigal
with `-p meta --min-gene 15 --min-edge-gene 15 --max-overlap 15`. The final
argument controls the number of assembly jobs running concurrently (default:
1); it is not passed as Pyrodigal's per-assembly thread setting.

Give every assembly a unique MAG or sample basename. For example, rename
sample-specific `final.contigs.fa` inputs to `<sample_id>.fa` before placing
them together. The script checks for conflicting output names and existing
protein outputs before starting prediction. Avoid whitespace and `|` in source
names because those characters complicate downstream sequence identifiers.

### 2. Filter peptides

```bash
python3 01_GFMSC/02_filter_smORFs.py \
  results/01_GFMSC/proteins/cFMD \
  results/01_GFMSC/filtered/cFMD_smORFs.faa

python3 01_GFMSC/02_filter_smORFs.py \
  results/01_GFMSC/proteins/grain \
  results/01_GFMSC/filtered/fermented_grain_smORFs.faa
```

Filtering is performed in this order:

1. Require `start_type=ATG` in the prediction header.
2. Remove terminal stop characters (`*`).
3. Remove one initial methionine (`M`), if present.
4. Retain processed peptide sequences of 5–50 amino acids, inclusive.

The length range therefore applies **after** initial methionine removal. The
script does not exclude partial ORFs or filter internal stop/ambiguous residues.
It reads protein files directly inside the specified directory, sorted by filename.

Optional arguments are `--min-length`, `--max-length` and `--allow-non-atg`.
The default example uses the intended ATG-only, processed 5–50-aa definition;
changing these arguments changes the input catalog.

Sequence IDs are written as `<source>|<original_protein_id>`. For example,
`sample_01.proteins.faa` containing `contig_12_3` becomes
`sample_01|contig_12_3`. The `.proteins.faa` suffix is removed. Ensure source
names are unique across both datasets and retain this convention in downstream
tables. It differs from older IDs that retained `.proteins` in the prefix.
Detailed prediction headers and coordinates are not retained in the filtered
FASTA; keep the original prediction files for traceability.

### 3. Build the combined catalog

```bash
bash 01_GFMSC/03_build_GFMSC.sh \
  results/01_GFMSC/filtered/cFMD_smORFs.faa \
  results/01_GFMSC/filtered/fermented_grain_smORFs.faa \
  results/01_GFMSC/catalog 64 64000
```

The script concatenates cFMD peptides followed by fermented-grain peptides and
clusters the combined collection once. The final two arguments specify the
number of threads and the memory limit in MB; defaults are 64 and 64,000.

| CD-HIT parameter | Setting |
| --- | --- |
| Sequence identity (`-c`) | 0.5 |
| Word length (`-n`) | 3 |
| Shorter-to-representative length ratio (`-s`) | 0.95 |
| Alignment coverage of the longer sequence (`-aL`) | 0.95 |
| Accurate assignment mode (`-g`) | 1 |
| Alignment information (`-p`) | 1 |
| Identifier truncation (`-d`) | 0 |

The `0.95` settings do not indicate 95% sequence identity. Preserve input order,
software versions and all clustering settings when comparing catalog versions.

| Output in `catalog/` | Description |
| --- | --- |
| `GFMSC_all_smORFs.faa` | Combined redundant peptide sequences |
| `GFMSC_representatives.faa` | One representative peptide per family |
| `GFMSC_representatives.faa.clstr` | Representative and member identifiers for each family |

The command preserves the recorded clustering parameters, including the absence
of an explicit short-sequence cutoff (`-l`). CD-HIT's documented default cutoff
is 10. Compare `.clstr` membership with the combined FASTA to determine whether
short inputs were excluded; do not assume every filtered peptide entered a family.
Family membership also does not mean that all members have identical sequences.

### 4. Assign contig taxonomy

```bash
bash 01_GFMSC/04_assign_contig_taxonomy.sh \
  data/fermented_grain_assemblies/sample_01.fa \
  databases/GTDB/gtdb \
  results/01_GFMSC/taxonomy/sample_01 \
  results/01_GFMSC/tmp/taxonomy_sample_01 64
```

Run this command for each fermented-grain assembly with a unique output prefix
and temporary directory. Use the same sample identifier as in prediction.
Keep assemblies separate if their contig IDs are reused across samples.

The script runs `mmseqs easy-taxonomy` and produces contig-level assignments.
To annotate smORFs, join each peptide to its source contig and transfer the
contig assignment. For cFMD-derived peptides, inherit taxonomy from the source
MAG metadata. **Neither transfer/join operation is implemented by this script.**
Retain source-aware contig IDs and original prediction records for those joins.

### 5. Search reference databases

```bash
bash 01_GFMSC/05_assess_smORF_novelty.sh \
  results/01_GFMSC/catalog/GFMSC_representatives.faa \
  results/01_GFMSC/novelty \
  results/01_GFMSC/tmp/novelty \
  databases/UniRef50/uniref50 \
  databases/SwissProt/swissprot \
  databases/NR/nr \
  databases/GMSC/gmsc
```

The script applies `mmseqs easy-search` to each supplied target using
`-e 0.001 --min-seq-id 0.3 --cov-mode 2 -c 0.8`, with additional recorded search
settings `--max-accept 1 --start-sens 1 --sens-steps 3 -s 7 --threads 64`.
Coverage is measured relative to the query. The thread count is fixed at 64 in
this script, rather than exposed as an argument.

Outputs are named `<database_basename>.m8`. Use distinct database basenames
(as in the example), since equal basenames would reuse output and temporary
paths. Keep the search targets consistent with the study database snapshots.

Each `.m8` is a tab-separated alignment table in the default BLAST-like format.
Column 1 contains the query ID. The script does not generate a final known/novel
table or calculate database overlap percentages. To obtain those summaries,
combine the unique query IDs with accepted hits across all four databases and
compare them with the full representative ID list. Peptides with no accepted hit
are described as lacking detectable homology **under the selected databases
and search thresholds**, rather than as absolutely novel sequences.

## Outputs for downstream analysis

Use `GFMSC_representatives.faa` as the input to bioactivity prediction in
`02_bioactive_peptides/`. Retain the combined FASTA, `.clstr`, source metadata
and original protein predictions to recover family membership and source records.

The manuscript reports 6,877,140 redundant smORFs and 2,237,774 families. These
are reference counts, not results automatically verified by these scripts.
Exact reproduction requires the original inputs, gene caller, filtering rules,
input order and database snapshots. The five scripts do not yet implement a
tabular representative/member/source mapping, taxonomy transfer or novelty
summary generation.

## Software documentation

- [Pyrodigal command-line interface](https://pyrodigal.readthedocs.io/en/v3.6.3/guide/cli.html)
- [CD-HIT user guide](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki)
- [MMseqs2 user guide](https://github.com/soedinglab/MMseqs2/wiki)
