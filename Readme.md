# Bulk_RNAseq_Pipeline_GSE158550

**Obtained raw data** from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305)  
**Reference paper** from [Nature communications](https://doi.org/10.1038/s41467-018-08133-6)

**Environment:** Windows Subsystem for Linux (Ubuntu)  
**System Requirements:** At least **16 GB RAM** (recommended)

---

## Project goal
This repository documents a **beginner-friendly, reproducible** Bulk RNA-Seq analysis workflow (download → QC → trimming → alignment → QC → counting → counts matrix preparation). The repository stores *commands exactly as executed* by the author, brief explanations, small result files (counts, QC HTMLs), scripts (text) and screenshots of terminal outputs so learners can reproduce the pipeline.

**Important:** Large raw files (FASTQ/SRA, BAMs, reference FASTA / indexes) are **not** committed to this repo. See **Data & reference policy** below.

---

## Data & reference policy (why we don’t commit everything)
- **Did not commit** raw `.sra`, raw `.fastq.gz`, `.bam` and aligner index files (they are huge).  
- Kept small, essential files in the repo: `sample_info.csv`, `featureCounts` counts matrix (CSV), a few representative QC HTMLs (MultiQC, FastQC) and small logs.  

---

## How to reproduce this analysis (high-level steps)
1. Clone repo and create environment.  
2. Download raw data & references by running the provided commands or scripts.  
3. Run QC → trimming → alignment → counting from the steps below.  
4. Copy the final counts matrix to `R_Pipeline/` and run downstream Differential Expression analysis using the 'DESeq2_tutorial_GSE106305.Rmd' file.  
5. Inspect example QC and Qualimap HTMLs included in `fastq/`, 'fastqc_results' & 'multiqc_report' and screenshots in `Screenshots/`.

---




## Step 1 - SRA fetch and FASTQ

```
#1
sudo apt install sra-toolkit
```

_Installs the SRA Toolkit on Debian/Ubuntu systems so you can use `prefetch`, `fastq-dump`, etc._

```
#2
prefetch SRR7179504
```

_Downloads the SRA file for the run `SRR7179504` into the local SRA cache._

```
#3
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ~/SRR7179504/SRR7179504.sra
```

_Converts the `.sra` file to gzipped FASTQ(s) and places them in the `fastq` folder; flags used: gzip output, skip technical reads, keep read IDs, keep only reads that pass filters, split paired reads correctly and apply SRA clipping._

```
#4
python3 fastq_download.py
```

_Runs a Python helper script (`fastq_download.py`) present in the root directory for downloading all the SRR files together and then converts '.sra' files to FASTQ gzips._
_Commands 2 and 3 can be skipped and this can be run directly for full automation._

```
#5
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz  > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz  > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz  > LNCAP_Hypoxia_S2.fastq.gz
```

_Concatenates four SRR FASTQ files to produce the combined sample files `LNCAP_Normoxia_S1.fastq.gz`, `LNCAP_Normoxia_S2.fastq.gz`, `LNCAP_Hypoxia_S1.fastq.gz` and `LNCAP_Hypoxia_S2.fastq.gz` respectively._

```
#6
mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz
```

_Renames `SRR7179536_pass.fastq.gz` to the biologically meaningful sample name `PC3_Normoxia_S1.fastq.gz`,`SRR7179537_pass.fastq.gz` to `PC3_Normoxia_S2.fastq.gz`, `SRR7179540_pass.fastq.gz` to `PC3_Hypoxia_S1.fastq.gz`and `SRR7179541_pass.fastq.gz` to `PC3_Hypoxia_S2.fastq.gz` respectively._

```
#7
rm SRR*
```

_Removes (deletes) all files starting with `SRR` in the current directory._  
**Note:** This permanently removes those files - keep a screenshot or backup if you want to keep the raw SRR files.




## Step 2 - FASTQC and MULTIQC

```
#8
sudo apt install fastqc
```

_Installs FastQC on Debian/Ubuntu systems so the `fastqc` program is available._

```
#9
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
```

_Runs FastQC on all FASTQ files in the folder `fastq/`, writes per-sample FastQC output into a new folder `fastqc_results/` and uses 8 threads to speed up processing._

```
#10
sudo apt install multiqc
```

_Installs MultiQC (system-level) so you can aggregate many QC reports into a single summary report._

```
#11
multiqc fastqc_results/ -o multiqc_report/
```

_Runs MultiQC over the `fastqc_results/` folder and writes a combined report into folder `multiqc_report/` (HTML and supporting files)._  
**MultiQC screenshot '11_multiqc.png' is available in '/Screenshots' folder for reference.**




## Step 3 - TRIMMOMATIC and FASTQC again

```
#12
java -jar trimmomatic-0.40.jar SE -threads 4 fastq/LNCAP_Hypoxia_S1.fastq.gz fastq/LNCAP_Hypoxia_S1_trimmed.fastq TRAILING:10 -phred33
```

**Download Trimmomatic from [here](http://www.usadellab.org/cms/?page=trimmomatic)**  
_Runs Trimmomatic in single-end (SE) mode on `fastq/LNCAP_Hypoxia_S1.fastq.gz`, writes trimmed output to `fastq/LNCAP_Hypoxia_S1_trimmed.fastq`, trims trailing bases with quality < 10, and assumes Phred+33 quality encoding._

```
#13
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
```

_Runs FastQC on all FASTQ files in `fastq/` (this will capture both original and trimmed files if both are present), writes FastQC reports to `fastqc_results/` and uses 8 threads._  
_We do this to check any differences between the original and trimmed sequences._  
**We only performed Trimming on one sample only just to check how much it affects the quality of the sequence, you may go ahead and perform QC on all the files.**




## Step 4 - Get Human Genome / Annotation

```
#14
wget ftp://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
```

_Downloads the Ensembl GTF annotation file (GRCh38 release 114) from the Ensembl FTP server._

```
#15
gunzip Homo_sapiens.GRCh38.114.gtf.gz
```

_Uncompresses the downloaded `.gtf.gz` to produce `Homo_sapiens.GRCh38.114.gtf`._




## Step 5 - Setting up and running Hisat2 & samtools using script

```
#16
sudo apt install hisat2
```

_Installs HISAT2 on Debian/Ubuntu so the `hisat2` binary is available._

```
#17
sudo apt install samtools
```

_Installs SAMtools (for SAM/BAM conversion, sorting, indexing and stats)._

```
#18
cd fastq
```

_Change directory into the `fastq/` folder where the sample FASTQ files and the alignment script 'hisat2alignment.sh' are located._

```
#19
./hisat2alignment.sh
```

_Runs the author's alignment script `hisat2alignment.sh`. This script should call `hisat2` for each sample and pipes output through `samtools` to create sorted/indexed BAM files and alignment summary statistics._  
**Script screenshot '19_hisat2.png' is available in '/Screenshots' folder for reference.**




## Step 6 - Install Qualimap & run RNA-Seq QC

Download Qualimap v2.3 from [here](http://qualimap.conesalab.org/) and keep in root directory.

```
#20
cd ..
```

_Change to the parent/root directory._

```
#21
sudo apt install r-base-core
```

_Installs R (core) which is needed by Qualimap dependencies and some R-based installation scripts._

```
#22
unzip qualimap_v2.3.zip
```

_Unzip Qualimap distribution into root directory._

```
#23
cd qualimap_v2.3
sudo Rscript scripts/installDependencies.r
cd ..
```

_Runs an R script to install any R dependencies required by Qualimap._

**Before going further, you must make changes to qualimap bash file to remove the default 1200M memory allocation, this will otherwise cause low memory problems in low end systems otherwise. In this case, we set the memory to 6G.**  
**You can check the error '23_error_qualimap.png' screenshot in '/Screenshots' folder for reference.**

```
#24
./qualimap_v2.3/qualimap rnaseq -bam fastq/LNCAP_Hypoxia_S1.bam -gtf Homo_sapiens.GRCh38.114.gtf –outdir fastq/gtf –java-mem-size=6G
./qualimap_v2.3/qualimap rnaseq -bam fastq/LNCAP_Hypoxia_S2.bam -gtf Homo_sapiens.GRCh38.114.gtf –outdir fastq/gtf –java-mem-size=6G
./qualimap_v2.3/qualimap rnaseq -bam fastq/LNCAP_Normoxia_S1.bam -gtf Homo_sapiens.GRCh38.114.gtf –outdir fastq/gtf –java-mem-size=6G
./qualimap_v2.3/qualimap rnaseq -bam fastq/LNCAP_Normoxia_S2.bam -gtf Homo_sapiens.GRCh38.114.gtf –outdir fastq/gtf –java-mem-size=6G
./qualimap_v2.3/qualimap rnaseq -bam fastq/PC3_Hypoxia_S1.bam -gtf Homo_sapiens.GRCh38.114.gtf –outdir fastq/gtf –java-mem-size=6G
./qualimap_v2.3/qualimap rnaseq -bam fastq/PC3_Hypoxia_S2.bam -gtf Homo_sapiens.GRCh38.114.gtf –outdir fastq/gtf –java-mem-size=6G
./qualimap_v2.3/qualimap rnaseq -bam fastq/PC3_Normoxia_S1.bam -gtf Homo_sapiens.GRCh38.114.gtf –outdir fastq/gtf –java-mem-size=6G
./qualimap_v2.3/qualimap rnaseq -bam fastq/PC3_Normoxia_S2.bam -gtf Homo_sapiens.GRCh38.114.gtf –outdir fastq/gtf –java-mem-size=6G
```

_Runs Qualimap RNA-seq on the `LNCAP_Hypoxia_S1.bam` using the provided '.gtf' and directs output to `fastq/gtf` with 6 GB Java heap size. The same Qualimap command is run for each sample BAM; see the script for all sample lines._
**You can visualize '.bam' files using [IGV Tools](https://igv.org/).**




## Step 7 - Feature counts with subread

```
#25
sudo apt install subread
```

_Installs the `subread` package that provides `featureCounts`._

```
#26
cd fastq
mkdir -p quants
```

_Creates an output directory `quants/` to store count results in the fastq directory._  

**Make sure to keep the Homo_sapiens.GRCh38.114.gtf file in fastq folder or define its path properly in script.**

```
#27
./featurecounts.sh
```

_Runs the `featurecounts.sh` script which calls `featureCounts` (subread) to generate the gene × sample counts matrix._  

**In case of permission denial error in #27, run the below command first.**

```
#27a
chmod +x featurecounts.sh
```
**You can see the error and execution '27_featurecounts.png' screenshot in 'Screenshots' folder for reference.**



## Step 8 - Generate counts matrix from FeatureCounts

```
#28
cd ..
sudo apt install -y python3-pandas
```

_Run this command if 'pandas' library is not available and return to the root directory._

```
#29
python3 countsmatrix_wholedata.py
```

_Runs the author's Python script which reads per-sample `featureCounts` outputs and writes a combined counts matrix 'GSE106305_counts_matrix.csv' in ~fastq/quants folder._




## Step 9 - Gene annotation for R & Cleaning counts

```
#30
cd fastq
wc -l Homo_sapiens.GRCh38.114.gtf
awk '$3 == "gene"' Homo_sapiens.GRCh38.114.gtf > genes_only.gtf
wc -l genes_only.gtf
```

_'wc' command counts the number of lines in the original '.gtf' (gives an idea of file size/complexity)._
_Main command here is awk, where we Filter the '.gtf' to keep only rows whose 3rd column is `gene`, writing the result to `genes_only.gtf`. This produces an annotation file with only gene-level records._

**After this step, we will open 'gene_only.gtf' in Excelsheet to arrange the data in 3 columns- Geneid, Genebiotype and Genesymbol. Save this file as 'GRCh38annotation.csv'.**




## Concepts involved and Decision-making

This section explains the reasoning and biological/computational concepts behind the main choices made in this pipeline. It’s meant to help learners understand **why** specific tools, parameters and workflows were chosen.

---

### Why use the latest genome annotation?
- The annotation file (e.g., *Homo_sapiens.GRCh38.114.gtf*) defines genomic features such as genes, exons and transcripts.
- Using the **latest Ensembl release** ensures:
  - Newly discovered genes and isoforms are included.
  - Deprecated or mis-annotated genes from older releases are corrected.
  - Gene biotypes and IDs are compatible with current databases (important for downstream tools like DESeq2 or Ensembl BioMart).
- Mismatched or outdated annotations can cause incomplete or incorrect gene counts.

---

### Why HISAT2 and not STAR?
- **HISAT2** was chosen for its:
  - **Low memory footprint** — suitable for mid-range systems (~16 GB RAM) and WSL setups.
  - **Fast and accurate alignment** of reads to the reference genome using FM-index and graph-based alignment for spliced reads.
- **STAR**, while faster for large servers, can require >32 GB RAM for human genomes, which is infeasible on many local systems.
- HISAT2 integrates seamlessly with **SAMtools** and **featureCounts** in lightweight pipelines.

---

### Why Trimming (Trimmomatic) before QC?
- Sequencing adapters and low-quality bases at the 3′ ends of reads can cause misalignments or poor mapping quality.
- **Trimmomatic** removes these, improving read quality.
- Running **FastQC** before and after trimming allows visual comparison of improvement in read quality metrics.

---

### Why use FastQC and MultiQC together?
- **FastQC** gives per-sample quality reports (per-base quality, GC content, adapter contamination).
- **MultiQC** aggregates all FastQC outputs into one HTML report for easy cross-sample comparison.

---

### Why perform QC again after alignment (Qualimap)?
- Pre-alignment QC tells us about raw data quality.
- **Qualimap RNA-seq** assesses *alignment quality*:
  - Proportion of reads mapping to exons, introns, intergenic regions.
  - Coverage uniformity across transcripts.
  - 3′/5′ bias and GC distribution.
- This helps detect potential library preparation issues or poor mapping.

---

### Why featureCounts instead of HTSeq-count?
- **featureCounts (Subread)** is faster, multi-threaded and memory-efficient.
- It can handle multiple '.bam' files in one run and uses simple, consistent column formatting.
- HTSeq-count is Python-based and single-threaded - better for small projects, but slower for large datasets.

---

### Why DESeq2 for Differential Expression?
- **DESeq2** models read counts with a **Negative Binomial Distribution**, accounting for biological variability and library size.
- It provides:
  - Normalization by size factors.
  - Shrinkage estimation for fold-change accuracy.
  - Integrated variance stabilization and visualization tools.
- Widely accepted as the gold standard for RNA-seq count-based analysis.

---

### Why not use TPM or FPKM for DE analysis?
- TPM/FPKM normalize for transcript length and sequencing depth but **do not model count dispersion**.
- DESeq2 works directly on raw counts - it performs internal normalization suitable for statistical testing.
- TPM/FPKM are better for expression *visualization*, not for differential testing.

---

### Why separate scripts for each step?
- Keeps the workflow **modular** - beginners can run or debug individual steps.
- Easier to replace one component (e.g., change Trimmomatic to fastp) without affecting others.
- Encourages understanding of each tool’s input/output.

---

### Why use a '.gtf' rather than a '.gff' file?
- GTF (Gene Transfer Format) is the preferred annotation format for tools like 'featureCounts' and 'DESeq2'.
- GFF3 is more general but less standardized in attribute columns; GTF provides consistent `gene_id` and `transcript_id` fields.

---

### Why merge SRR runs per condition before alignment?
- GEO studies often provide multiple technical replicates per condition (e.g., SRR7179504–SRR7179507).
- Merging them creates a single, larger FASTQ for that condition → more reads, better coverage.
- Biological replicates (e.g., LNCAP_Hypoxia_S1 vs S2) are still kept separate for DESeq2 statistical testing.

---

### Why remove intermediate files (e.g., SRR*, .fastq.gz)?
- To save disk space and keep the working directory clean.
- Once concatenated and renamed, the intermediate 'SRR' files are no longer needed.

---

### Why not include raw data in GitHub?
- FASTQ, BAM and index files are extremely large (>10GB).
- GitHub’s file size limit is 100 MB per file.
- Instead, links and commands for downloading are provided in `Readme.md`.

---

### Why use Ensembl rather than RefSeq annotations?
- Ensembl provides well-maintained release cycles and consistent naming (gene IDs compatible with Bioconductor packages).
- RefSeq is equally reliable but sometimes lags behind in version updates.
- For reproducibility, always document which source and release number were used.

---

### Why re-check normalization and PCA before DE?
- PCA and sample clustering validate if replicates group together correctly.
- If a sample clusters far from others of the same condition, it might be an outlier (batch effect or technical error).

---

### Key Computational Concepts Highlighted
- **Alignment vs. pseudoalignment:** HISAT2 performs alignment; tools like Salmon/Kallisto do pseudoalignment (faster but count-level only).  
- **Count-based modeling:** DESeq2 uses raw integer counts, not normalized values.  
- **Quality Control feedback loop:** Each QC step informs decisions about trimming, filtering or re-mapping.  
- **Annotation consistency:** Using the same '.gtf' for featureCounts and downstream DE ensures gene IDs align perfectly.

---

### Summary
This pipeline prioritizes **clarity, reproducibility and accessibility**:
- Each step is transparent (manual commands instead of black-box scripts).
- Choices are guided by system limitations and reproducible best practices.
- Beginners can swap components as they learn (e.g., try STAR or Salmon later) and compare results.

---