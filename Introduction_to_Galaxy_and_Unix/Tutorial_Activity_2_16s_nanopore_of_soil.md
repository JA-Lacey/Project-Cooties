# Soil 16S Nanopore Analysis Tutorial
# Jake A Lacey

**This tutorial is modified from the galaxy project tutorials to provide a unix/command-line version of the tutorial.** 


https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/nanopore-16S-metagenomics/tutorial.html#re-evaluate-datasets-quality


This tutorial guides you through a complete workflow for analysing soil-derived 16S rRNA Nanopore sequencing data. You will learn how to:

1. **Set up a reproducible environment** with all necessary bioinformatics tools using Conda.
2. **Organise your directory structure** for raw and processed data, QC outputs, and taxonomy results.
3. **Download Nanopore reads** for different soil compartments.
4. **Generate a manifest** file for batch processing.
5. **Perform quality control (QC)** on raw and processed reads using FastQC and MultiQC.
6. **Trim and filter reads** using Porechop and Fastp.
7. **Compress processed files** to save space.
8. **Classify reads taxonomically** with Kraken2 against a SILVA-based database.
9. **Visualise taxonomic profiles** interactively with Krona.

---

## 1. Environment Setup

Create a Conda environment to isolate dependencies and make the workflow reproducible.

```bash
# Create Conda environment with required tools
conda create -n Soil16S_env \
  fastqc multiqc fastp porechop kraken2 krona taxonkit gnu-parallel -y

# Verify and activate
conda info --envs
conda activate Soil16S_env

```
**This tutorial is modified from the galaxy project tutorials to provide a unix/command-line version of the tutorial.** 

https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/beer-data-analysis/tutorial.html#solution-8 
```

```
**This tutorial is modified from the galaxy project tutorials to provide a unix/command-line version of the tutorial.** 

https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/beer-data-analysis/tutorial.html#solution-8 
```

*Tools installed:*

- **fastqc**: raw read QC
- **multiqc**: aggregates QC reports
- **porechop**: Nanopore adapter trimming
- **fastp**: quality filtering and length trimming
- **kraken2**: taxonomic classification
- **krona**: interactive plots
- **taxonkit**: taxonomy lineage formatting
- **gnu-parallel**: parallel execution of commands

---

## 2. Directory Structure

Organise your project under a root directory (e.g. `/Soil_16S_Nanopore`). Create subdirectories for raw reads, QC, processed reads, taxonomy, and visualisation:

```bash
bash
CopyEdit
# Define project root
export PROJECT_ROOT=/Soil_16S_Nanopore

# Create hierarchy
mkdir -p $PROJECT_ROOT/data/reads/
mkdir -p $PROJECT_ROOT/analysis/qc/fastqc_raw/
mkdir -p $PROJECT_ROOT/analysis/qc/porechop/
mkdir -p $PROJECT_ROOT/analysis/qc/fastp/
mkdir -p $PROJECT_ROOT/analysis/qc/fastqc_processed/
mkdir -p $PROJECT_ROOT/data/reads_processed/
mkdir -p $PROJECT_ROOT/analysis/taxonomy/kraken2/
mkdir -p $PROJECT_ROOT/analysis/visualisation/krona

# Navigate to reads folder
cd $PROJECT_ROOT/data/reads/

```

Resulting structure:

```
arduino
CopyEdit
/Soil_16S_Nanopore
├── data
│   ├── reads           # raw downloads
│   └── reads_processed # filtered/trimmed files
└── analysis
    ├── qc
    │   ├── fastqc_raw
    │   ├── porechop
    │   ├── fastp
    │   └── fastqc_processed
    ├── taxonomy
    │   └── kraken2
    └── visualisation
        └── krona

```

---

## 3. Downloading Nanopore Reads

Retrieve four soil samples (bulk and rhizosphere, top and bottom horizons):

```bash
bash
CopyEdit
wget https://zenodo.org/record/4274812/files/bulk_bottom.fastq.gz
wget https://zenodo.org/record/4274812/files/bulk_top.fastq.gz
wget https://zenodo.org/record/4274812/files/rhizosphere_bottom.fastq.gz
wget https://zenodo.org/record/4274812/files/rhizosphere_top.fastq.gz

```

Move or symlink these into your raw reads directory:

```bash
bash
CopyEdit
mv *.fastq.gz $PROJECT_ROOT/data/reads/
cd $PROJECT_ROOT/data/reads/

```

---

## 4. Manifest File Creation

Generate a tab-delimited manifest (`reads_paths.tab`) mapping each sample ID to its file path:

```bash
bash
CopyEdit
# List all FASTQ files with full path
find $PROJECT_ROOT/data/reads/ -name "*.fastq.gz" -exec realpath {} \; > reads_paths.txt

# Extract sample ID (e.g. bulk_bottom) and pair with path
awk -F '/' '{ sub(/\.fastq\.gz$/, "", $NF); print $NF "\t" $0 }' reads_paths.txt \
  > reads_paths.tab

```

Inspect with:

```bash
bash
CopyEdit
head reads_paths.tab
# Expected format:
# bulk_bottom   /Soil_16S_Nanopore/data/reads/bulk_bottom.fastq.gz
# ... etc.

```

---

## 5. QC of Raw Reads

### 5.1 FastQC in a Loop

Run FastQC on each sample, outputting to `analysis/qc/fastqc_raw/`:

```bash
bash
CopyEdit
while read sample path; do
  fastqc "$path" -o $PROJECT_ROOT/analysis/qc/fastqc_raw/"$sample" -t 8
done < reads_paths.tab

```

### 5.2 Aggregate with MultiQC

Combine all FastQC reports into one HTML:

```bash
bash
CopyEdit
cd $PROJECT_ROOT/analysis/qc/fastqc_raw/
multiqc . -o ../multiqc_raw

```

---

## 6. Trimming and Filtering

### 6.1 Adapter Removal with Porechop

```bash
bash
CopyEdit
# In parallel
cat $PROJECT_ROOT/data/reads/reads_paths.tab \
  | parallel -j 1 --colsep '\t' \
      'porechop -i {2} -o $PROJECT_ROOT/analysis/qc/porechop/{1}_trim.fastq -t 8'

```

Or equivalently:

```bash
bash
CopyEdit
while read sample path; do
  porechop -i "$path" \
           -o $PROJECT_ROOT/analysis/qc/porechop/"${sample}_trim.fastq" \
           -t 8
done < reads_paths.tab

```

### 6.2 Quality Filtering with Fastp

```bash
bash
CopyEdit
while read sample path; do
  fastp \
    -i $PROJECT_ROOT/analysis/qc/porechop/"${sample}_trim.fastq" \
    --disable_adapter_trimming \
    --qualified_quality_phred 9 \
    --length_required 1000 \
    --length_limit 2000 \
    --disable_trim_poly_g \
    -o $PROJECT_ROOT/analysis/qc/fastp/"${sample}_filt.fastq"
done < reads_paths.tab

```

---

## 7. QC of Processed Reads

1. **Build new manifest** for filtered reads:
    
    ```bash
    bash
    CopyEdit
    realpath $PROJECT_ROOT/analysis/qc/fastp/*.fastq > processed_reads_paths.txt
    awk -F '/' '{ sub(/\.fastq$/, "", $NF); print $NF "\t" $0 }' \
        processed_reads_paths.txt > processed_reads.tab
    
    ```
    
2. **FastQC on filtered reads**:
    
    ```bash
    bash
    CopyEdit
    mkdir -p $PROJECT_ROOT/analysis/qc/fastqc_processed/
    while read sample path; do
      fastqc "$path" -o $PROJECT_ROOT/analysis/qc/fastqc_processed/"$sample" -t 8
    done < processed_reads.tab
    
    ```
    
3. **MultiQC aggregation**:
    
    ```bash
    bash
    CopyEdit
    cd $PROJECT_ROOT/analysis/qc/fastqc_processed/
    multiqc . -o ../multiqc_processed
    
    ```
    

---

## 8. Compressing Processed Reads

Navigate to the processed reads folder and gzip files to save space:

```bash
bash
CopyEdit
cd $PROJECT_ROOT/analysis/qc/fastp/
gzip *.fastq

```

*Optionally* move compressed files into your long-term storage directory:

```bash
bash
CopyEdit
mv *.fastq.gz $PROJECT_ROOT/data/reads_processed/

```

---

## 9. Taxonomic Classification with Kraken2

Use a SILVA-based Kraken2 database located at `/home/mdu/resources/kraken2/rdp`.

```bash
bash
CopyEdit
cat processed_reads.tab \
  | parallel -j 4 --colsep '\t' \
      'kraken2 {2} \
         --db /home/mdu/resources/kraken2/rdp \
         --threads 8 \
         --use-names \
         --confidence 0.1 \
         --report $PROJECT_ROOT/analysis/taxonomy/kraken2/{1}_report.txt \
         --output $PROJECT_ROOT/analysis/taxonomy/kraken2/{1}_output.txt'

```

---

## 10. Preparing Krona Inputs

Convert Kraken2 outputs into Krona-compatible counts:

```bash
bash
CopyEdit
# Extract taxids and counts
awk '$1 == "C" { for (i=1;i<=NF;i++) if ($i~/taxid/) { gsub(/[()]/,"",$(i+1)); print $(i+1); break } }' \
    $PROJECT_ROOT/analysis/taxonomy/kraken2/*_output.txt \
  | sort | uniq -c > taxid_counts.txt

# Generate lineage TSV
awk '{print $2}' taxid_counts.txt \
  | taxonkit lineage \
  | taxonkit reformat -f "{k};{p};{c};{o};{f};{g};{s}" \
  > lineage.tsv

# Create Krona input table
paste <(awk '{print $1}' taxid_counts.txt) \
      <(cut -f2 lineage.tsv | tr ';' '\t') \
  > krona_input.txt

```

Finally, build the Krona chart:

```bash
bash
CopyEdit
ktImportText krona_input.txt \
  -o $PROJECT_ROOT/analysis/visualisation/krona/soil16S_krona.html

```

Open `soil16S_krona.html` in your browser to explore the taxonomic composition interactively.

