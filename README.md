# DBiT-seq_FFPE
This repository aims to share the raw data processing and visualization codes used in DBiT-seq_FFPE project.

<img src="https://github.com/edicliuyang/Spatial_FFPE/blob/master/scheme.png">

**Refer to our lastest preprint paper for more details: https://www.biorxiv.org/content/10.1101/2020.10.13.338475v1.abstract**

## Abstract

DBiT-seq was a spatial multiomic sequencing technique that can also be directly applied for FFPE samples. We have successfully demonstrated spatial transcriptome sequencing of embryonic and adult mouse FFPE tissue sections at cellular level (25μm pixel size) with high coverage (>1,000 genes per pixel). Spatial transcriptome of a E10.5 mouse embryo identified all major anatomical features in the brain and abdominal region. 

### Raw Fastq data processing

We did the illumina pair-end 100 sequencing using Hiseq 4000. We pool two samples for each sequencing lane. 

**The DBiT-seq Raw fastq file:**

<p><img src="https://github.com/MingyuYang-Yale/DBiT-seq/blob/master/Pre-processing/schematic.png" alt="foo bar" title="train &amp; tracks" /></p>

Read 1: Contains the RNA sequence

Read 2: Contains the sptial Barcode and UMIs



For data preprocessing and generation of expression matrix from raw data, please refer to https://github.com/MingyuYang-Yale/DBiT-seq

This repository includes the main R scripts used for the visualization of the sequencing data, including clustering, Differential expression gene analysis and cell annotation with reference scRNA-seq data. 

One MATLAB script to automatic identify pixels on tissue was also included.

