# MOSAIT: Multi-Omics Spatial Analysis Integration Tool

## 1. Functional Overview

MOSAIT is a unified command-line toolkit for spatial transcriptomics data processing.
It integrates multiple pipelines and provides standardized preprocessing, analysis, and visualization.

The toolkit consists of four modules:

- run: Native pipeline for flexible raw data processing
- zumi: Wrapper for zUMIs pipeline
- astro: Wrapper for ASTRO pipeline
- visual: 

Core functionalities:

- Demultiplexing (barcode / UMI extraction)
- Integration with external pipelines (zUMIs / ASTRO)
- Standardized workflow execution
- 

---

## 2. Installation Guide

### 2.1 install by conda

### 2.2 clone

### 2.3 install zUMIs and Astro

---

## 3. Modules Description

### 3.1 run (Core Pipeline)

Custom pipeline for processing raw FASTQ data.

```bash
tool run --config config.yaml
```

### 3.2 zUMIs (zUMIs wrapper)

Wrapper for the zUMIs pipeline.

```bash
tool zumi --config config.yaml
```

### 3.3 ASTRO (ASTRO wrapper)

Wrapper for the zUMIs pipeline.

```bash
tool astro --config config.yaml
```

### 3.4 visual (Visualization)

Generate spatial spot maps on microscopy images.

```bash
yourtool visual \
    --in tissue.png \
    --out position.txt \
    --p 50 \
    --t 0
```

## 4.Parameter Description

MOSAIT

## Parameters

| Parameter Name | Required | Relevant Steps | Default | Description |
|----------------|----------|----------------|---------|-------------|
| Parameter Name | Required | Relevant Steps | Default | Description |
