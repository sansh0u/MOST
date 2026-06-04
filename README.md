# MOST: Multi-Omics Spatial Tool

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

---
## 2. Getting start

### 2.1 Installation
MOST is also available on bioconda. Thus you can easily install MOST with Conda.
```bash
conda install -c conda-forge most
```
 Check if the installation was successful. If you see the help documentation, the installation is complete:
```bash
most --help
```


### 2.2 Install zUMIs and Astro

ASTRO and zUMIs can be obtained either by cloning their GitHub repositories or by downloading the packaged versions from the links provided below.

### ASTRO
- GitHub: https://github.com/gersteinlab/ASTRO#
- Direct download: https://
```bush
```
### zUMIs
- GitHub: https://github.com/sdparekh/zUMIs
- Direct download: https://
```bush
```

---

## 3. Modules Description

### 3.1 `run` (Native Pipeline)

The ``` run ``` module is the native workflow implemented in MOST for processing raw sequencing data.

Depending on the value of the ```method``` parameter in the configuration file, the workflow automatically selects the appropriate processing pipeline:

| Method | Supported Technologies | Backend  
|----------------|----------|----------------|
| RNA | DBiT-seq, Patho-DBiT | ST-Pipeline | 
| ATAC| Spatial ATAC, Spatial CUT&Tag, Patho-ATAC, Patho-CUT&Tag  | Chromap | 

| Method | Technology | Backend |
|----------|----------|----------|
| RNA | DBiT-seq | ST-Pipeline |
| RNA | Patho-DBiT | ST-Pipeline |
| ATAC | Spatial ATAC | Chromap |
| ATAC | Spatial CUT&Tag | Chromap |
| ATAC | Patho-ATAC | Chromap |
| ATAC | Patho-CUT&Tag | Chromap |

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
| -------------- | -------- | ------- | ----------------------------------------------- |
| `project`      | Yes      | `Project       | Project name used for output files and logs.    |
| `method`       | Yes      | -   | Analysis mode. Supported values: `RNA`, `ATAC`. |
| `out_dir`      | Yes      | -       | Output directory for all generated results.     |
| `threads`      | No       | `16`     | Number of CPU threads used during analysis.     |

### Input Files

| Parameter Name        | Required | Default | Description                                                                  |
| --------------------- | -------- | ------- | ---------------------------------------------------------------------------- |
| `sequence_file.file1` | Yes      | -       | Read 1 FASTQ file, typically containing biological sequences for alignment.             |
| `sequence_file.file2` | Yes      | -       | Read 2 FASTQ file, typically containing barcode and UMI sequences. |

### Reference Files

| Parameter Name           | Required | Default       | Description                                                                                              |
| ------------------------ | -------- | ------------- | -------------------------------------------------------------------------------------------------------- |
| `reference.index_file`   | Yes      | -             | Genome index used by Chromap alignment steps.                                                              |
| `reference.fa_file`      | Yes      | -             | Reference genome FASTA file.                                                                             |
| `reference.gtf_file`     | Yes      | -             | Gene annotation GTF file.                                                                                |
| `reference.star_index`   | Yes      | -             | STAR genome index directory.                                                                             |
| `reference.barcode_file` | No       | [Default barcode whitelist](most/barcode/20240614_2500barcode_AB_update.txt) | Spatial barcode whitelist file. If not specified, MOST uses the bundled barcode whitelist distributed with the package (most/barcode/). |

### Advanced Parameters

These parameters are optional and should only be modified when processing custom library designs.

| Parameter Name              | Required | Default          | Description                                                                        |
| --------------------------- | -------- | ---------------- | ---------------------------------------------------------------------------------- |
| `advanced.primer`           | No       | `CAAGCGTTGGCTTCTCGCATCT` | Primer sequence used for library construction.                                     |
| `advanced.linker1`          | No       | `GTGGCCGATGTTTCGCATCGGCGTACGACT` | Linker 1 sequence.                                                                 |
| `advanced.linker2`          | No       | `ATCCACGTGCTTGAGAGGCCAGAGCATTCG` | Linker 2 sequence.                                                                 |
| `advanced.UMI`              | No       | Auto-detected  | Position of the UMI sequence in Read 1. Coordinates are 1-based and inclusive.     |
| `advanced.BC`               | No       | Auto-detected  | Position of spatial barcodes in Read 1. Multiple barcode regions can be specified. |
| `advanced.hdist`            | No       | `3`              | Maximum Hamming distance allowed when matching barcodes.                           |
| `advanced.adapter`          | No       | `AAGCAGTGGTATCAACGCAGAGTGAATGGG`    | Adapter sequence to trim from reads.                                               |
| `advanced.adapter_mismatch` | No       | `1`              | Maximum number of mismatches allowed during adapter detection.                     |
---

### Notes

* Positions specified in `UMI` and `BC` use 1-based coordinates.

  Example:

  ```yaml
  UMI: (23-32)
  BC: (33-40,71-78)
  ```

