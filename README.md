# Transcriptomics of Wood Frog and American Bullfrog *Bd* Response

This repository contains code, data, and figures that support:

Eskew, E.A., B.C. Shock, E.E.B. LaDouceur, K. Keel, M.R. Miller, J.E. Foley, and B.D. Todd. 2018. [Gene expression differs in susceptible and resistant amphibians exposed to *Batrachochytrium dendrobatidis*.](https://doi.org/10.1098/rsos.170910) Royal Society Open Science 5: 170910.

The primary aim of this project was to help elucidate the mechanisms by which some amphibian host species are resistant to or tolerant of the fungal pathogen *Batrachochytrium dendrobatidis* (*Bd*) given that infection in many species results in high mortality, leading to population declines and species-level extinctions. Thus, we experimentally exposed two common amphibian species expected to differ in disease response, wood frogs and American bullfrogs, to *Bd* and measured a suite of outcomes. This repository provides data and code related to frog survival, body mass change, and *Bd* infection prevalence and load (evaluated via both qPCR and histology). Ultimately, these data support the key product from our work: an RNA sequencing dataset generated from ventral skin tissues from 87 study animals. RNA-seq samples were collected from both amphibian host species at multiple time points following *Bd* exposure. After sequencing and read quality processing, reads were mapped to [an existing amphibian reference transcriptome](http://datadryad.org/resource/doi:10.5061/dryad.j6676) using [Salmon](http://salmon.readthedocs.io/en/latest/). This repository includes Salmon data resulting from two alternative read processing pipelines: 1) processing using a combination of [ConDeTri](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0026314) and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) or 2) processing via the [expHTS pipeline](https://github.com/msettles/expHTS). Ultimately, we use the latter data for downstream analyses. Differential gene expression analyses were conducted using the [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) package.

Raw sequence data for this project are available at NCBI's Sequence Read Archive accession number [SRP111327](https://www.ncbi.nlm.nih.gov/sra/?term=SRP111327).  

--- 

### Repository Structure

- `data\` contains all data files necessary to conduct the analyses
	- `annotation\` subdirectory contains gene annotation data used in `bioinformatics_analyses.R`
	- `expHTS_preprocess\` subdirectory contains metadata resulting from expHTS processing of raw RNA-seq reads
	- `salmon\` subdirectory contains read count data (generated from the Salmon software) used in `bioinformatics_analyses.R`
- `figures\` contains all figures that are output from analysis scripts
- `R\` contains code for functions that are sourced in the analysis scripts
- `scripts\` contains all scripts necessary to recreate this project's analyses
	- `bioinformatics_analyses.R` contains code for RNA-seq data analyses
	- `body_measurement_analyses.R` contains code for frog body mass and length analyses
	- `expHTS_preprocess_analyses.R` contains code for plotting of expHTS metadata and calculation of post-processing RNA-seq read counts
	- `histology_analyses.R` contains code for analyses comparing infection quantification with histology and qPCR
	- `infection_analyses.R` contains code for infection prevalence and load analyses (using qPCR data)
	- `survival_analyses.R` contains code for frog survival analyses

---

### Reproducing the Analyses

Assuming you have all the required R packages, recreating analyses and figures is as simple as running the relevant script from the `scripts` directory. Scripts are written such that figures will be output to the `figures` directory. Note that some scripts, particuarly `bioinformatics_analyses.R`, will take quite a while to work through.

---

### License

This project is licensed under the terms of the MIT license.
