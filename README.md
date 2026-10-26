# Environmental Resistome and Homogenization of Extremotolerant Communities in Jet Propulsion Laboratory Cleanroom Environments

[![NASA OSDR OSD-944](https://img.shields.io/badge/NASA%20OSDR-OSD--944-blue.svg)](https://osdr.nasa.gov/bio/repo/data/studies/OSD-944)
[![DOI](https://img.shields.io/badge/DOI-10.26030%2F3zww--vy60-green.svg)](https://doi.org/10.26030/3zww-vy60)
[![License: CC0-1.0](https://img.shields.io/badge/License-CC0--1.0-lightgrey.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-brightgreen.svg)](https://www.python.org/)

**Author:** Md Abu Sadat Soash  
**Affiliation:** Department of Genetic Engineering and Biotechnology, Shahjalal University of Science and Technology, Sylhet, Bangladesh  
**Repository:** [https://github.com/Soash/nasa-osd944-resistome-analysis](https://github.com/Soash/nasa-osd944-resistome-analysis)

---

## Abstract

Jet Propulsion Laboratory (JPL) cleanroom environments adhere to strict planetary protection protocols to minimize microbial bioburden, inadvertently creating selective pressures that drive the survival of extremotolerant, multidrug-resistant (MDR) microbes. This study maps the spatial distribution and ecological diversity of the environmental resistome within a JPL Mars Exploration Program cleanroom using short-read metagenomic sequencing data (**OSD-944**). Following the removal of mock community controls and normalization of sequencing depth to Coverage Per Million (CPM), the resistome was characterized against the KEGG and BacMet databases. 

Spatial analysis revealed the highest localized burdens of antimicrobial resistance (AMR) in human transition zones, specifically the floor by cabinets (1,211.6 CPM) and gowning rooms (1,066.5 CPM), while critical spacecraft assembly surfaces exhibited near-zero AMR loads. Despite these observational differences, statistical evaluations including Kruskal-Wallis ($H = 8.5686, p = 0.285$) and PERMANOVA (pseudo-$F = 1.123, p = 0.321$) indicated no significant spatial variance in total AMR burden or microbial community composition across the facility footprint. Taxonomic profiling identified primary carriers of resistance mechanisms as *Acinetobacter*, *Brevundimonas*, and *Paracoccus* species, alongside eukaryotic fungi such as *Aspergillus*. The detection of broad-spectrum efflux pumps, such as *mexJK*, and dual-resistance genes highlights a uniform, highly selected microbiome shaped by stringent sterilization protocols.

---

## Key Findings & Metrics

### 1. Assembly Success & Master Functional Matrix
From 38 raw samples, 32 samples processed under the standardized "Prep 3" Illumina Nextera XT protocol were curated. Environmental cleanroom samples demonstrated a 100% assembly rate, yielding a unified matrix of 62,195 functional gene annotations.

| Metric | Environmental Samples | Control Samples | Total Dataset |
| :--- | :---: | :---: | :---: |
| **Targeted Samples (Prep 3)** | 18 | 14 | 32 |
| **Failed Assembly** | 0 | 7 | 7 |
| **Successfully Assembled** | 18 | 7 | 25 |
| **Assembly Success Rate** | **100%** | 50.0% | **78.1%** |
| **Total Functional Genes (ORFs)** | 56,282 | 5,913 | 62,195 |
| **AMR Genes Detected (KEGG)** | 447 | 62 | 509 |
| **Dual-Resistance Genes (AMR + BRG)** | **64** | 0 | **64** |

### 2. Cleanroom Resistome & Dual-Resistance Determinants
* **Total ARGs & BRGs:** 509 Antimicrobial Resistance Genes (KEGG ko01504) and 904 Biocide/Metal Resistance Genes (BacMet) were mapped across the dataset. Isolating the cleanroom samples confirmed 447 bona fide environmental ARGs.
* **Dual Resistance (Cross-Resistance):** 64 genes confer simultaneous resistance to both biocides/metals and clinical antibiotics (e.g., *merR* family [K19591], *tetR*/*acrR* family transcriptional repressors [K18301], and *mexJK* broad-spectrum efflux systems). All 64 dual-resistance genes were localized exclusively within environmental samples.
* **Primary Microbial Reservoirs:**
  * ***Acinetobacter* (863.1 CPM):** Mediated predominantly by streptothricin acetyltransferases and TetR/AcrR family repressors.
  * ***Brevundimonas* (424.5 CPM):** Enriched in fluoroquinolone resistance proteins, BlaI family penicillinase repressors, and PhoQ sensor histidine kinases.
  * ***Paracoccus* (407.8 CPM):** Characterized by aminocyclitol acetyltransferases and multidrug efflux pumps.
  * ***Aspergillus* (20.1 CPM):** Eukaryotic fungal lineage carrying resistance factors (preliminary detection via prokaryotic ORF calling).

### 3. Spatial Distribution & Facility Homogenization
* **Localized Burdens:** Human transition zones demonstrated the highest AMR densities:
  * Floor by cabinets: **1,211.6 CPM**
  * Gowning room floor: **1,066.5 CPM**
  * Spacecraft hardware assembly surfaces: **Near-zero CPM**
* **Statistical Homogeneity:**
  * **Kruskal-Wallis $H$-test:** $H = 8.5686, p = 0.285$ (no significant spatial difference in total AMR burden).
  * **PERMANOVA (Beta Diversity, 8 locations):** pseudo-$F = 1.123, p = 0.321$ (no significant divergence in microbial community structure).
  * **Conclusion:** Cleanroom decontamination regimens act as an overarching selective filter that homogenizes the resident extremotolerant community across the entire facility footprint.

---

## Dataset Information

* **Repository:** NASA Open Science Data Repository (OSDR)
* **Study ID:** [OSD-944](https://osdr.nasa.gov/bio/repo/data/studies/OSD-944) (GeneLab ID: `GLDS-775`)
* **Study Title:** *JPL Mars Exploration Program Interlab Cleanroom Metagenomics*
* **DOI:** [10.26030/3zww-vy60](https://doi.org/10.26030/3zww-vy60)
* **Sequencing Platform:** Illumina NovaSeq X (paired-end 2×150 bp, Nextera XT library prep, Prep 3)
* **Upstream Processing:** NASA OSDR pipeline (FastQC, fastp, Bowtie2 decontamination, Kraken2 human-read filtering, MEGAHIT assembly, Prodigal ORF calling, KOFamScan, CAT/BAT taxonomic classification).

---

## Pipeline Architecture

The end-to-end reproducible analysis pipeline is implemented across 6 sequential Jupyter Notebooks:

```text
OSD-944 ISA Metadata & S3 Annotations
                   │
                   ▼
       1_clean_metadata.ipynb          ──► Filter Prep 3 swabs vs. controls
                   │
                   ▼
       2_clean_assaydata.ipynb         ──► Process sequencing depths & metrics
                   │
                   ▼
     3_download_annotations.ipynb      ──► Fetch gene coverage & taxonomic TSVs
                   │
                   ▼
     4_build_master_matrix.ipynb       ──► Merge & CPM-normalize (62,195 ORFs)
                   │
                   ▼
           5_resistome.ipynb           ──► Map KEGG AMR, BacMet BRGs, & Dual Genes
                   │
                   ▼
             6_stat.ipynb              ──► Alpha/Beta diversity & PERMANOVA tests
```

### Notebook Descriptions

1. **`1_clean_metadata.ipynb`**: Parses raw ISA-Tab archives (`s_OSD-944.txt`, `i_Investigation.txt`). Filters for short-read low-biomass Illumina sequencing (`Prep 3`) and isolates environmental cleanroom swabs from controls. Outputs `data/clean_metadata.tsv`.
2. **`2_clean_assaydata.ipynb`**: Extracts sequencing metrics including raw read counts, post-filtering depth, DNA concentration, and human/control decontamination rates. Outputs `data/clean_assaydata.tsv`.
3. **`3_download_annotations.ipynb`**: Automated batch downloader fetching `gene-coverage-annotation-and-tax.tsv` files directly from NASA OSDR AWS S3 buckets for all target samples into `data/annotations/`.
4. **`4_build_master_matrix.ipynb`**: Consolidates 25 individual assembly tables, applies Coverage Per Million (CPM) normalization relative to sample sequencing depths, and builds `data/master_functional_matrix.tsv`.
5. **`5_resistome.ipynb`**: Cross-references functional genes against the KEGG Antimicrobial Resistance database (`ko01504`) and the BacMet database. Identifies antibiotic resistance classes, biocide/metal resistance, and 64 dual-resistance determinants. Exports `data/env_amr_df.tsv` and generates Figures 1–3.
6. **`6_stat.ipynb`**: Computes Alpha diversity (Shannon index) and Beta diversity (Bray-Curtis dissimilarity, PCoA). Runs non-parametric Kruskal-Wallis and PERMANOVA tests across 8 cleanroom spatial locations. Generates Figures 4–5.

---

## Repository Structure

```text
.
├── 1_clean_metadata.ipynb           # Step 1: Parse ISA metadata and sample filters
├── 2_clean_assaydata.ipynb          # Step 2: Extract sequencing metrics & depth
├── 3_download_annotations.ipynb     # Step 3: Automated S3 download of OSDR tables
├── 4_build_master_matrix.ipynb      # Step 4: Consolidate annotations & CPM normalization
├── 5_resistome.ipynb                # Step 5: Profiling ARGs, BRGs, & dual resistance
├── 6_stat.ipynb                     # Step 6: Statistical testing & diversity metrics
├── data/
│   ├── clean_metadata.tsv           # Filtered sample metadata
│   ├── clean_assaydata.tsv          # Curated sequencing metrics
│   ├── master_functional_matrix.tsv # CPM-normalized 62,195-gene master matrix
│   ├── env_amr_df.tsv               # Environmental AMR profile dataframe
│   ├── database/                    # KEGG and BacMet reference mapping files
│   ├── annotations/                 # Downloaded sample annotation TSVs
│   └── result/                      # Publication-ready figures
│       ├── amr_by_location.png      # Figure 1: AMR burden across cleanroom locations
│       ├── amr_by_genus.png         # Figure 2: AMR contribution by microbial genus
│       ├── acinetobacter_mechanisms.png # Figure 3: Functional resistance in Acinetobacter
│       ├── alpha_diversity.png      # Figure 4: Alpha diversity by sampling location
│       └── beta_diversity_true_pcoa.png # Figure 5: PCoA ordination of community structures
├── nasa-osd944-resistome-analysis.docx # Full research manuscript draft
├── supplementary.txt                # Study protocols, parameters, and metadata notes
└── README.md
```

---

## Installation & Environment Setup

### 1. Clone Repository
```bash
git clone https://github.com/Soash/nasa-osd944-resistome-analysis.git
cd nasa-osd944-resistome-analysis
```

### 2. Create and Activate Virtual Environment

**Windows (PowerShell):**
```powershell
python -m venv venv
.\venv\Scripts\Activate.ps1
```

**macOS / Linux:**
```bash
python3 -m venv venv
source venv/bin/activate
```

### 3. Install Dependencies
```bash
pip install --upgrade pip
pip install pandas requests matplotlib seaborn scipy scikit-bio
```

---

## Reproduction Guide

Execute the Jupyter Notebooks in numbered numerical order:

```bash
jupyter notebook
```
1. Run **`1_clean_metadata.ipynb`** to generate `data/clean_metadata.tsv`.
2. Run **`2_clean_assaydata.ipynb`** to generate `data/clean_assaydata.tsv`.
3. Run **`3_download_annotations.ipynb`** to pull raw OSDR tables from AWS S3.
4. Run **`4_build_master_matrix.ipynb`** to assemble `data/master_functional_matrix.tsv`.
5. Run **`5_resistome.ipynb`** to perform resistome mapping and export resistance figures.
6. Run **`6_stat.ipynb`** to run statistical hypothesis tests (Kruskal-Wallis, PERMANOVA) and diversity plots.

---

## References

1. **Guan, L., Hinzer, A., Parker, C., & Cooper, M. (2026).** *JPL Mars Exploration Program Interlab Cleanroom Metagenomics*, NASA Open Science Data Repository, Version 1, [DOI: 10.26030/3zww-vy60](http://doi.org/10.26030/3zww-vy60).
2. **Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M., & Ishiguro-Watanabe, M. (2023).** KEGG for taxonomy-based analysis of pathways and genomes. *Nucleic Acids Research*, 51(D1), D587–D592.
3. **Pal, C., Bengtsson-Palme, J., Rensing, C., Kristiansson, E., & Larsson, D. G. J. (2014).** BacMet: antibacterial biocide and metal resistance genes database. *Nucleic Acids Research*, 42(D1), D737–D743.
4. **La Duc, M. T., Dekas, A., Osman, S., Moissl, C., Newcombe, D., & Venkateswaran, K. (2007).** Isolation and characterization of bacteria capable of tolerating the extreme conditions of clean room environments. *Applied and Environmental Microbiology*, 73(8), 2600–2611.
5. **The scikit-bio development team. (2020).** *scikit-bio: A bioinformatics library for Python*.

---

## License & Acknowledgements

* **Data License:** NASA Open Science Data Repository (OSDR) data is provided under the [Creative Commons CC0 1.0 Universal Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/).
* **Code License:** Open for academic and research use under CC-BY (Attribution).
* **Acknowledgements:** The author acknowledges the National Aeronautics and Space Administration (NASA) and the GeneLab / Open Science Data Repository team for public access to the metagenomic sequencing archives.
