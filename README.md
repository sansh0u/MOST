# MOST: Multi-Omics Spatial Tool
## Table of Contents
- [1. Functional Overview](#1-functional-overview)
- [2. Getting Started](#2-getting-started)
- [3. Modules Description](#3-modules-description)
- [4. Parameter Description](#4-parameter-description)



## 1. Functional Overview

MOST is a unified command-line toolkit for spatial transcriptomics data processing.
It integrates multiple preprocessing pipelines into a standardized framework and simplifies spatial omics data analysis for users with different levels of computational experience.

Key features:
- Flexible raw FASTQ processing workflow
- Automatic barcode and UMI position detection
- Integration with ASTRO and zUMIs pipelines
- Standardized preprocessing and workflow execution
- Visualization utilities for spatial data
- Built-in default parameters to reduce manual configuration

![Workflow Diagram](/workflow.png)

## 2. Getting start

### 2.1 Installation
MOST is available on conda-forge. Thus you can easily install MOST with Conda.
```bash
conda install -c conda-forge most
```
 Check if the installation was successful. If you see the help documentation, the installation is complete:
```bash
most --help
```


### 2.2 Install zUMIs and Astro

Installation of zUMIs and ASTRO is not required for MOST. However, MOST includes integrated interfaces to these tools, enabling users to perform RNA expression quantification with zUMIs and Patho-DBiT analysis with ASTRO directly from MOST. 

ASTRO and zUMIs can be obtained either by cloning their GitHub repositories or by downloading the packaged versions from the links provided below. 

### ASTRO
- GitHub: https://github.com/gersteinlab/ASTRO#
- Direct download: https://
```bash
#clone the repository
git clone git@github.com:gersteinlab/ASTRO.git 
#enter the directory named "python"
cd python
#install dependencies and build/install
pip install .
#check if the installation was successful
ASTRO --help
```


### zUMIs
- GitHub: https://github.com/sdparekh/zUMIs
- Direct download: https://
```bash
#clone the repository
git clone https://github.com/sdparekh/zUMIs.git
```

### 2.3 Download Reference Files
MOST requires appropriate reference files for RNA and ATAC data processing. Users can download the reference files provided by the MOST project or prepare their own custom reference files.

Human (GRCh38 / hg38)

```bash
mkdir reference/hg38
cd reference/hg38
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz
```
Mouse (GRCm38 / mm10)
```bash
mkdir reference/mm10
cd reference/mm10
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
```
### 2.4 Create Reference Index

2.4.1 Chromap Index

Chromap requires a genome index before alignment. The index only needs to be generated once for each reference genome.
```bash
# Example for human hg38
cd reference/hg38
#decompress the reference genome while keeping the original .gz file
gunzip -c GRCm38.primary_assembly.genome.fa.gz > GRCm38.primary_assembly.genome.fa
#build Chromap index
chromap -i -r GRCm38.primary_assembly.genome.fa -o chromap_index
```
2.4.2 STAR Index

STAR requires a separate genome index for RNA-seq alignment.
```bash
# Example for mouse mm10
cd reference/mouse
gunzip -c GRCm38.primary_assembly.genome.fa.gz > GRCm38.primary_assembly.genome.fa

gunzip -c gencode.vM23.annotation.gtf.gz > gencode.vM23.annotation.gtf

STAR --runMode genomeGenerate \
    # STAR genome index directory
    --genomeDir star_index \
    --genomeFastaFiles GRCm38.primary_assembly.genome.fa \
    --sjdbGTFfile gencode.vM25.annotation.gtf \
    --sjdbOverhang 99 \
    --runThreadN 16
```
## 3. Modules Description

### 3.1 `run` (Native Pipeline)

The ``` run ``` module is the native workflow implemented in MOST for processing raw sequencing data.

Depending on the value of the ```method``` parameter in the configuration file, the workflow automatically selects the appropriate processing pipeline:


| Method | Technology | Backend |
|:----------|:--------:|:-------|
| RNA | DBiT-seq | ST-Pipeline |
| RNA | Patho-DBiT | ST-Pipeline |
| ATAC | Spatial ATAC | Chromap |
| ATAC | Spatial CUT&Tag | Chromap |
| ATAC | epi-Patho-ATAC | Chromap |
| ATAC | epi-Patho-CUT&Tag | Chromap |

Example:

```bash
most run --config config.yaml
```

### 3.2  `zumis` (zUMIs Wrapper)

Wrapper for the zUMIs pipeline. Designed for RNA spatial transcriptomics data, it provides a fast and standardized interface to run zUMIs within MOST.

Unlike the native zUMIs workflow, the zumis module uses the MOST configuration file (```config.yaml```) and automatically generates the required zUMIs settings.

Before using this module, zUMIs must be installed separately. During the first run, provide the path to the `zUMIs.sh` script:

```bash
most zumis --config config.yaml --l /path/to/zUMIs.sh
```

The zUMIs path will be saved in the MOST configuration and does not need to be specified again in subsequent runs.

After initialization, the workflow can be executed using the standard MOST configuration file:

```bash
most zumis --config config.yaml
```
### 3.3 `astro` (ASTRO wrapper)

Wrapper for the ASTRO pipeline. ASTRO is a flexible workflow for whole-transcriptome spatial omics data analysis, supporting both coding and non-coding RNA profiling. Originally developed for FFPE-based Patho-DBiT datasets, it has also been validated on multiple spatial transcriptomics platforms. Through MOST, users can run ASTRO with simplified configuration and automated parameter handling.



```bash
most astro --config config.yaml
```

### 3.4 `visual` (Visualization)

Generate spatial spot maps on microscopy images.

```bash
most visual \
    --in tissue.png \
    --out position.txt \
    --p 50 \
    --t 0
```


## 4. Parameter Description

The MOST workflow is configured through a unified YAML file that is shared across most modules. The following table describes all supported parameters.

### General Parameters

| Parameter Name | Required | Default | Description                                     |
|:----------|:--------:|:-------:|:------------|
| `project`      | Yes      | Project       | Project name used for output files and logs.    |
| `method`       | Yes      | -   | Analysis mode. Supported values: `RNA`, `ATAC`. |
| `out_dir`      | Yes      | -       | Output directory for all generated results.     |
| `threads`      | No       | `16`     | Number of CPU threads used during analysis.     |

### Input Files

| Parameter Name        | Required | Default | Description                                                                  |
|:----------|:--------:|:-------:|:------------|
| `sequence_file.file1` | Yes      | -       | Read 1 FASTQ file, typically containing biological sequences for alignment.             |
| `sequence_file.file2` | Yes      | -       | Read 2 FASTQ file, typically containing barcode and UMI sequences. |

### Reference Files

| Parameter Name           | Required | Default       | Description                                                                                              |
|:----------|:--------:|:-------:|:------------|
| `reference.index_file`   | Yes      | -             | Genome index used by Chromap alignment steps.                                                              |
| `reference.fa_file`      | Yes      | -             | Reference genome FASTA file (e.g., genome.fa.gz).                                                                             |
| `reference.genome`      | Yes      | -             | Reference genome assembly name used for QC plot generation (e.g., hg38).                                                                |
| `reference.gtf_file`     | Yes      | -             | Gene annotation GTF file (e.g., genome.gtf.gz).                                                                                |
| `reference.star_index`   | Yes      | -             | STAR genome index directory.                                                                             |
| `reference.barcode_file` | No       | [Default barcode whitelist](most/barcode/20240614_2500barcode_AB_update.txt) | Spatial barcode whitelist file. If not specified, MOST uses the bundled barcode whitelist distributed with the package (most/barcode/). |

### Advanced Parameters

These parameters are optional and should only be modified when processing custom library designs.

| Parameter Name  | Required | Default | Description |
|:----------|:--------:|:-------:|:------------|
| `advanced.primer`           | No       | `CAAGCGTTGGCTTCTCGCATCT` | Primer sequence used for library construction.                                     |
| `advanced.linker1`          | No       | `GTGGCCGATGTTTCGCATCGGCGTACGACT` | Linker 1 sequence.                                                                 |
| `advanced.linker2`          | No       | `ATCCACGTGCTTGAGAGGCCAGAGCATTCG` | Linker 2 sequence.                                                                 |
| `advanced.UMI` | No       | Auto-detected  | Position of the UMI sequence in Read 1. Coordinates are 1-based and inclusive.     |
| `advanced.BC`               | No       | Auto-detected  | Position of spatial barcodes in Read 1. Multiple barcode regions can be specified. |
| `advanced.hdist` | No       | `3`              | Maximum Hamming distance allowed when matching barcodes.                           |
| `advanced.adapter`          | No       | `AAGCAGTGGTATCAACGCAGAGTGAATGGG`    | Adapter sequence to trim from reads.                                               |
| `advanced.adapter_mismatch` | No       | `1`              | Maximum number of mismatches allowed during adapter detection.                     |


### Notes

* Positions specified in `UMI` and `BC` use 1-based coordinates.

  Example:

  ```yaml
  UMI: (23-32)
  BC: (33-40,71-78)
  ```

