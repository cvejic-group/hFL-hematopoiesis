# Disomy Hi-C analysis

These scripts process the two disomy libraries, **16216 and 17759**, and produce the interaction calls, QC tables and figures.

Run them in this order:

1. **01_run_hicup.py** - align reads to hg38, filter invalid contacts and remove duplicates.
2. **02_make_contact_maps.py** - convert each library's BAM file into a Juicer `.hic` contact map.
3. **03_call_interactions.R** - call interactions with HiC-DC+ at 5, 10 and 25 kb; save call counts and distance plots.
4. **04_library_qc.R** - report valid pairs, duplicates, and short-range cis, long-range cis and trans contacts.
5. **05_replicate_qc.R** - compare the two libraries at 25 kb: correlations, significant-contact overlap, Jaccard index and Fisher tests.

## Before running

Use an environment with Python 3, HiCUP, Bowtie2, Samtools, Juicer Tools/Java, and R with HiCDCPlus, BSgenome.Hsapiens.UCSC.hg38, data.table, ggplot2 and scales.

Copy `config.example.json` to `config.json` and enter the output directory, hg38 FASTA path and resource settings. Copy `samples.example.tsv` to `samples.tsv` and enter the actual paired FASTQ paths. Use one row per lane; lanes from the same library share the same donor number. Keep pilot and study libraries separate.

Activate the environment before running. On Iris, this is usually:

```bash
conda activate hicup_clean
```

## Run everything

Place this folder on Iris at:

`mardera1@islogin01.mskcc.org:/data1/offitk/mardera1/analysis/disomy_lihic_scripts/`

Then submit:

```bash
sbatch /data1/offitk/mardera1/analysis/disomy_lihic_scripts/run_all.sh \
  /data1/offitk/mardera1/analysis/disomy_lihic_scripts \
  /data1/offitk/mardera1/analysis/disomy_lihic_scripts/config.json \
  /data1/offitk/mardera1/analysis/disomy_lihic_scripts/samples.tsv
```

`run_all.sh` runs steps 01-05 and saves software versions. It uses an explicit code path so Slurm can find the scripts. Outputs go into the directory specified in `config.json`. Use a new output directory for a fresh run.

Each script also gives its command-line usage when run without arguments. Existing `.hic` maps can be used directly in step 03; existing 25-kb interaction tables can be used directly in step 05.

## Reading the results

Interaction significance uses **q <= 0.05** with HiC-DC+'s within-chromosome BH adjustment. Pearson correlations use log(1 + counts); Spearman correlations use observed counts. The scripts explain these choices in comments.

The replicate-QC calculations were checked against the original disomy tables. The full raw-data workflow still needs to run in the cluster environment. Historical HiCUP/HiC-DC+ versions and the previously quoted 219M/189M totals still require the original run records.
