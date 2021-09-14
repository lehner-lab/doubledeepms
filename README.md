Welcome to the GitHub repository for the following publication: [Global mapping of the energetic and allosteric landscapes of protein binding domains (Faure AJ, Domingo J & Schmiedel JM et al., 2021)]()

Here you'll find an R package with all scripts to reproduce the figures and results from the computational analyses described in the paper.

# Table Of Contents

* **1. [Required Software](#required-software)**
* **2. [Installation Instructions](#installation-instructions)**
* **3. [Required Data](#required-data)**
* **4. [Pipeline Modes](#pipeline-modes)**
* **5. [Pipeline Stages](#pipeline-stages)**

# Required Software

To run the doubledeepms pipeline you will need the following software and associated packages:

* **[_R_](https://www.r-project.org/) >=v3.6.1** (bio3d, Biostrings, coin, Cairo, data.table, ggplot2, GGally, hexbin, plot3D, reshape2, RColorBrewer, ROCR, stringr, ggrepel)

The following software is optional:

* **[_Python_](https://www.python.org/) v3.8.6** (pandas, numpy, matplotlib, tensorflow, scikit-learn) 
* **[DiMSum](https://github.com/lehner-lab/DiMSum) v1.2.8** (pipeline for pre-processing deep mutational scanning data i.e. FASTQ to counts)

# Installation Instructions

Open R and enter:

```
# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("lehner-lab/doubledeepms")

# Load
library(doubledeepms)

# Help
?doubledeepms
```

# Required Data

Fitness scores, thermodynamic models, pre-processed data and required miscellaneous files should be downloaded from [here](https://www.dropbox.com/s/w45fuv36o0uxe7p/Data.zip?dl=0) and unzipped in your project directory (see '_base_dir_' option) i.e. where output files should be written.

# Pipeline Modes

There are a number of options available for running the doubledeepms pipeline depending on user requirements.

* ## Basic (default)

Default pipeline functionality ('_startStage_' = 1) uses prefit thermodynamic models and fitness scores from DMS experiments (already processed with MoCHI and DiMSum respectively; see [Required Data](#required-data)) to reproduce all figures in the publication.

* ## Thermodynamic model inference with MoCHI

Pipeline stage 0 ('doubledeepms_fit_thermo_model') fits thermodynamic models to DMS data for the specified domains ('_tmodel_protein_'), using all available data or subsets of phenotypes/variants ('_tmodel_subset_'). Parallel computing using job arrays is reccommended while running monte carlo simluations to determine confidence intervals of model-inferred free energies ('_tmodel_job_number_'). **Note:** this stage can be resource intensive (up to 48h with 30GB of RAM for GB1). 

* ## Raw read processing

Raw read processing is not handled by the doubledeepms pipeline. FastQ files ([GSE184042](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184042)) from paired-end sequencing of replicate deep mutational scanning (DMS) libraries before ('input') and after selection ('output') were processed using **[DiMSum](https://github.com/lehner-lab/DiMSum)** (Faure and Schmiedel et al., 2020).

DiMSum command-line arguments and Experimental design files required to obtain variant counts from FastQ files are available [here](https://www.dropbox.com/sh/pfn0bbk14k27zvj/AABQnRykP9vtIhQFoXK_cWhJa?dl=0).

# Pipeline Stages

The top-level function **doubledeepms()** is the recommended entry point to the pipeline and by default reproduces the figures and results from the computational analyses described in the following publication: [Global mapping of the energetic and allosteric landscapes of protein binding domains (Faure AJ, Domingo J & Schmiedel JM et al., 2021)](). See [Required Data](#required-data) for instructions on how to obtain all required data and miscellaneous files before running the pipeline.

## Stage 0: Fit thermodynamic models with MoCHI

This stage ('doubledeepms_fit_thermo_model') fits thermodynamic models to variant fitness data from (_ddPCA_) DMS.

## Stage 1: Evaluate thermodynamic model results

This stage ('doubledeepms_thermo_model_results') evaluates thermodynamic model results and performance including comparing to literature _in vitro_ measurements (related to Figure 2).

## Stage 2: Add 3D structure metrics

This stage ('doubledeepms_structure_metrics') annotates single mutant inferred free energies with PDB structure-derived metrics.

## Stage 3: Fitness plots

This stage ('doubledeepms_fitness_plots') plots fitness distributions and scatterplots (related to Figure 1).

## Stage 4: Plot fitness heatmaps

This stage ('doubledeepms_fitness_heatmaps') plots single mutant fitness heatmaps (related to Figure 1).

## Stage 5: Plot free energy scatterplots

This stage ('doubledeepms_free_energy_scatterplots') plots single mutant free energy scatterplots (related to Figure 3).

## Stage 6: Plot free energy heatmaps

This stage ('doubledeepms_free_energy_heatmaps') plots single mutant free energy heatmaps (related to Figure 3).

## Stage 7: Protein stability plots

This stage ('doubledeepms_protein_stability_plots') produces protein stability plots (related to Figure 4).

## Stage 8: Binding interface plots

This stage ('doubledeepms_interface_mechanisms') produces binding free energy heatmaps for selected GRB2-SH3 residues (related to Figure 5).

## Stage 9: Allostery plots

This stage ('doubledeepms_allostery_plots') produces allostery plots (related to Figure 5 and Figure 6).

## Stage 10: Allostery scatterplots

This stage ('doubledeepms_allostery_scatterplots') produces free energy scatterplots of major allosteric sites and mutations (related to Figure 5 and Figure 6).

## Stage 11: Downsampling analysis

This stage ('doubledeepms_downsampling_analysis') evaluates thermodynamic model results and performance after downsampling (related to Figure 2).


