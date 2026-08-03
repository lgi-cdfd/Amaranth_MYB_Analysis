# Comprehensive Genome-Wide Curation, Evolutionary Characterization, and Transcriptional Dynamics of MYB Transcription Factors in *Amaranthus hypochondriacus*

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21770784-blue.svg)](https://doi.org/10.5281/zenodo.21770784)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

## Overview
This repository contains the custom computational pipelines, Python scripts, curated sequence datasets, and differential expression analysis workflows supporting the manuscript:

> **"Genome-wide characterization of MYB transcription factors in *Amaranthus hypochondriacus* uncovers evolutionary conservation and stress-dependent transcriptional regulation"**  
>  
> **Ramesh Eerapagula¹**, **Rakesh Singh²**, **Ajay Kumar Mahato¹\***  
>  
> ¹ Laboratory of Genome Informatics, BRIC-Centre for DNA Fingerprinting and Diagnostics (CDFD), Hyderabad, Telangana, India.  
> ² ICAR-National Bureau of Plant Genetic Resources, Pusa Campus, New Delhi, India.  
>  
> \* **Corresponding Author:** Ajay Kumar Mahato (`akmahato@cdfd.org.in`)

---

## Repository Structure

* `scripts/`: Custom Python and R scripts used for MYB superfamily curation, repeat counting, domain validation, DESeq2 differential expression, and Ka/Ks calculation.
* `data/`: Curated protein FASTA, coding sequence (CDS) FASTA, gene metadata, and RNA-seq count matrices across 24 biological libraries.
* `tables/`: Supplementary dataset tables containing physicochemical properties, MEME motif annotations, PlantCARE promoter cis-elements, and DEG statistics under abiotic (salinity, drought) and biotic (bacterial infection, insect herbivory) stresses.

---

## Curation Workflow & Methodology

The MYB transcription factor superfamily in *Amaranthus hypochondriacus* (Genome v2.0, Phytozome v13) was identified and classified using a multi-stage curation pipeline:

1. **Candidate Identification:** Cross-referencing HMMER v3.3.2 searches against Pfam-A (PF00249, PF08914), InterProScan v5.68-100.0, and PlantTFDB master records.
2. **Isoform Selection & Non-TF Filtering:** Retaining only the longest representative isoform per locus and systematically excluding SANT domain chromatin remodellers (PF01382) and response regulators (PF00072).
3. **Repeat Counting & Domain Merging:** Single-repeat Pfam HMM profiles (PF00249, PF08914, PF12776, PF13837, PF13873) were used for repeat counting. Overlapping domain intervals (minimum 15 residues) were merged to prevent duplicate counting.
4. **Structural Validation:** Tryptophan (W) residue conservation was manually verified across HTH repeats to confirm functional DNA-binding capacity.
5. **Subfamily Classification:** Classifying proteins into 1R-MYB (81), R2R3-MYB (43), 3R-MYB (3), and 5R-MYB (1) subfamilies.

---

## Usage Instructions

### Prerequisites
* Python 3.8 or higher
* `biopython >= 1.79`
* `pandas >= 1.3.0`
* `numpy >= 1.21.0`
* `R >= 4.2.0` (`DESeq2`, `pheatmap`, `ggplot2`)

### Running Curation Scripts
```bash
# Clone the repository
git clone https://github.com/lgi-cdfd/Amaranth_MYB_Analysis.git
cd Amaranth_MYB_Analysis

# Execute candidate integration and repeat counting
python3 scripts/01_candidate_myb_integration.py
python3 scripts/02_isoform_selection_and_repeat_counting.py
```

### Running Differential Expression Analysis
```bash
# Execute DESeq2 analysis across developmental and stress libraries
Rscript scripts/05_dge_star_deseq2_analysis.R
```

---

## Data Availability
* **Reference Genome & Baseline Transcriptomes:** NCBI SRA BioProject [PRJNA263128](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA263128)
* **Stress RNA-Seq Libraries:** NCBI SRA BioProject [PRJNA65409](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA65409)
* **Zenodo Archive DOI:** [10.5281/zenodo.21770784](https://doi.org/10.5281/zenodo.21770784)

---

## Citation
If you use the datasets or scripts in this repository, please cite:
> Eerapagula R, Singh R, Mahato AK. (2026). Genome-wide characterization of MYB transcription factors in *Amaranthus hypochondriacus* uncovers evolutionary conservation and stress-dependent transcriptional regulation. *BMC Genomics* (In Press).

---

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.