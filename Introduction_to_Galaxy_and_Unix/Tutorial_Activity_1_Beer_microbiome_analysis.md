# Beer Microbiome Analysis Tutorial
# Author Jake A Lacey 

**This tutorial is modified from the galaxy project tutorials to provide a unix/command-line version of the tutorial.** 


https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/beer-data-analysis/tutorial.html#solution-8 

This tutorial walks you through the complete workflow for analysing beer microbiome sequencing data. You will learn how to:

1. **Set up a reproducible environment** with all necessary bioinformatics tools using Conda.
2. **Organize the directory structure** for raw and processed data.
3. **Download sequencing reads** and create manifest files for batch processing.
4. **Perform quality control** on raw reads with FastQC and MultiQC.
5. **Trim and filter reads** using Porechop and Fastp.
6. **Generate QC reports** for processed reads.
7. **Run taxonomic classification** with Kraken2 and summarise results.
8. **Visualise taxonomic profiles** using Krona.

---

## 1. Setting up the Required Environments/Tools

To ensure consistency and ease of installation, we will use a Conda environment. This isolates dependencies and allows easy sharing of the setup.

For more information on install conda (miniconda) see link 
/docs/getting-started/miniconda/install#using-miniconda-in-a-commercial-setting 

```
# Create a new Conda environment named "Beer_microbiome_env" with all required tools
conda create -n Beer_microbiome_env multiqc fastqc fastp porechop kraken2 krona taxonkit gnu-parallel -y
```

- **multiqc**: aggregates reports from multiple QC tools into a single HTML report.
- **fastqc**: performs quality checks on raw sequencing reads.
- **fastp**: trims and filters reads based on quality metrics.
- **porechop**: adapter trimming specifically for Oxford Nanopore reads.
- **kraken2**: taxonomic classification using k-mer matching.
- **krona**: interactive visualisation of taxonomic data.
- **taxonkit**: efficient manipulation of NCBI taxonomy.
- **gnu-parallel**: parallelises shell commands across multiple cores.

```
# Verify the environment was created
conda info --envs

# Activate the environment
conda activate Beer_microbiome_env
```

After activation, you can confirm each tool is available by running, for example, `fastqc -h` and `porechop --help`.

---

## 2. Setting up the Directory Structure on the Server

Organise your project into separate folders for raw data, QC outputs, and final analysis. This promotes clarity and reproducibility.

```
# Create base directories for raw reads and QC analysis
mkdir -p /Beer_Microbiome/data/reads/
mkdir -p /Beer_Microbiome/analysis/qc/
```

it will look something like this is visual format

```kotlin
/Beer_Microbiome
├── data
│   └── reads
└── analysis
    └── qc
```

Move into the reads directory where we will download our data:

```
cd /Beer_Microbiome/data/reads/
```

---

## 3. Downloading Reads and Creating Manifest Files

1. **Download** the raw fastq files using `wget`:
    
    ```
    wget https://zenodo.org/record/7093173/files/ABJ044_c38189e89895cdde6770a18635db438c8a00641b.fastq
    ```
    
2. **Generate a manifest** of all `.fastq` files with their absolute paths:
    
    ```
    find . -type f -name "*.fastq" -exec realpath {} \; > reads_paths.txt
    ```
    
3. **Format the manifest** into a tab-delimited file (`reads_paths.tab`) with two columns: sample ID and file path. We extract the sample ID from the filename prefix (e.g., `ABJ044`):
    
    ```
    awk -F'/' '{ match($NF, /^([A-Z0-9]+)_/, m); print m[1] "\t" $0 }' reads_paths.txt > reads_paths.tab
    ```
    
4. **Inspect the manifest** using `less` or `head`:
    
    ```
    less reads_paths.tab
    ```
    
5. **Optionally** move intermediate files to a scratch directory for cleanup:
    
    ```
    mkdir -p /Beer_Microbiome/scratch/
    mv reads_paths.txt /Beer_Microbiome/scratch/
    ```
    

---

## 4. Quality Control of Raw Reads with FastQC and MultiQC

Perform individual QC with FastQC, then aggregate results with MultiQC.

```
# Run FastQC on each sample in parallel (1 at a time here, adjust -j for concurrency)
cat reads_paths.tab \
  | parallel -j 1 --colsep '\t' 'fastqc {2} -o /Beer_Microbiome/analysis/qc/fastqc_raw -t 8'

# Aggregate all FastQC outputs into a single report
tcd /Beer_Microbiome/analysis/qc/fastqc_raw
multiqc . -o ../multiqc_raw
```

- `o` specifies output directory.
- `t 8` uses 8 threads per FastQC job (adjust to match your system).

---

## 5. Improving Reads by Trimming and Filtering

### 5.1 Adapter and Primer Removal with Porechop

```
# Create output folder for trimmed reads
mkdir -p /Beer_Microbiome/analysis/qc/porechop/

# Run Porechop in parallel
cat reads_paths.tab \
  | parallel -j 1 --colsep '\t' 'porechop -i {2} -o /Beer_Microbiome/analysis/qc/porechop/{1}_trim.fastq -t 8'
```

### 5.2 Quality Filtering with Fastp

```
# Create output folder for filtered reads
mkdir -p /Beer_Microbiome/analysis/qc/fastp/

# Example for a single sample; loop or parallelise as needed
fastp \
  -i /Beer_Microbiome/analysis/qc/porechop/ABJ044_trim.fastq \
  --disable_adapter_trimming \
  --qualified_quality_phred 10 \
  --disable_trim_poly_g \
  -o /Beer_Microbiome/analysis/qc/fastp/ABJ044_trim_filt.fastq
```

- `-qualified_quality_phred 10` discards bases with quality < Q10.
- `-disable_trim_poly_g` prevents trimming poly-G tails common in Illumina but not desired here.

---

## 6. QC of Processed Reads

Re-run FastQC on the trimmed/filtered reads and aggregate with MultiQC:

```
mkdir -p /Beer_Microbiome/analysis/qc/fastqc_processed/

# FastQC on processed reads
fastqc /Beer_Microbiome/analysis/qc/fastp/ABJ044_trim_filt.fastq \
  -o /Beer_Microbiome/analysis/qc/fastqc_processed -t 8

# Aggregate with MultiQC
tcd /Beer_Microbiome/analysis/qc/fastqc_processed
multiqc . -o ../multiqc_processed
```

---

## 7. Organising Processed Reads and New Manifest

1. **Copy processed reads** back to the data directory for long-term storage:
    
    ```
    mkdir -p /Beer_Microbiome/data/reads_processed/
    cp /Beer_Microbiome/analysis/qc/fastp/ABJ044_trim_filt.fastq /Beer_Microbiome/data/reads_processed/
    ```
    
2. **Generate a new manifest** for processed reads:
    
    ```
    cd /Beer_Microbiome/data/reads_processed/
    find . -type f -name "*filt.fastq" -exec realpath {} \; > reads_trim_paths.txt
    awk -F'/' '{ match($NF, /^([A-Z0-9]+)_/, m); print m[1] "\t" $0 }' reads_trim_paths.txt > reads_trim_paths.tab
    ```
    

---

## 8. Taxonomic Classification with Kraken2

1. **Prepare output directories**:
    
    ```
    mkdir -p /Beer_Microbiome/analysis/taxonomy/kraken2/
    ```
    
2. **Run Kraken2** to classify reads and generate reports:
    
    ```
    # Classification report only
    cat reads_trim_paths.tab \
      | parallel -j 1 --colsep '\t' \
        'kraken2 {2} \
         --use-names \
         --threads 8 \
         --report /Beer_Microbiome/analysis/taxonomy/kraken2/{1}_report.txt'
    
    # Full output plus report
    cat reads_trim_paths.tab \
      | parallel -j 1 --colsep '\t' \
        'kraken2 {2} \
         --use-names \
         --threads 8 \
         --report /Beer_Microbiome/analysis/taxonomy/kraken2/{1}_report.txt \
         --output /Beer_Microbiome/analysis/taxonomy/kraken2/{1}_output.txt'
    ```
    
3. **Filter the report** to remove low-count or unclassified lineages:
    
    ```
    # Keep taxa with >5 reads
    awk '$2 > 5' \
        /Beer_Microbiome/analysis/taxonomy/kraken2/ABJ044_report.txt \
      > ABJ044_report_filtered.txt
    
    # Remove unclassified entries as well
    awk '$2 > 5 && $6 != "unclassified"' \
        /Beer_Microbiome/analysis/taxonomy/kraken2/ABJ044_report.txt \
      > ABJ044_report_classified.txt
    ```
    

---

## 9. Visualisation with Krona

Convert the Kraken2 output into Krona-compatible format and generate an interactive HTML:

```
# 1. Extract and count taxids from classified reads
awk '$1 == "C" { for (i=1; i<=NF; i++) if ($i ~ /taxid/) { gsub(/[()]/, "", $(i+1)); print $(i+1); break } }' \
    /Beer_Microbiome/analysis/taxonomy/kraken2/ABJ044_output.txt \
  | sort | uniq -c > taxid_counts.txt

# 2. Retrieve full taxonomic lineage and reformat
awk '{print $2}' taxid_counts.txt \
  | taxonkit lineage \
  | taxonkit reformat -f "{k};{p};{c};{o};{f};{g};{s}" \
  > lineage.tsv

# 3. Prepare Krona input (count + each taxonomic rank as a separate column)
paste <(awk '{print $1}' taxid_counts.txt) \
      <(cut -f2 lineage.tsv | tr ';' '\t') \
  > krona_input.txt

# 4. Generate the Krona chart
ktImportText krona_input.txt -o ABJ044_krona.html
```

- Open `ABJ044_krona.html` in a web browser to explore the taxonomic composition interactively.