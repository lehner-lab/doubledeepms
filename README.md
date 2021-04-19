Welcome to the GitHub repository for the following publication: "A comprehensive map of allosteric sites in a protein domain"

Here you'll find an R package with all scripts to reproduce the figures and results from the computational analyses described in the paper.

# Table Of Contents

* **1. [Required Software](#required-software)**
* **2. [Installation Instructions](#installation-instructions)**
* **3. [Required Data](#required-data)**
* **4. [Pipeline Modes](#pipeline-modes)**
* **5. [Pipeline Stages](#pipeline-stages)**

# Required Software

To run the doubledeepms pipeline you will need the following software and associated packages:

* **[_R_](https://www.r-project.org/) >=v3.6.1** (data.table, ggplot2, hexbin, optparse, parallel, reshape2, RColorBrewer)

The following packages are optional:

* **[_DiMSum_](https://github.com/lehner-lab/DiMSum)** 
* **[_MoCHI_](https://github.com/lehner-lab/MoCHI)** 

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

Fitness scores, thermodynamic models, pre-processed data and required miscellaneous files should be downloaded from [here]() to your project directory (see '_base_dir_' option) i.e. where output files should be written, and unzipped.

# Pipeline Modes

There are a number of options available for running the dimsumms pipeline depending on user requirements.

* ## Basic (default)

Default pipeline functionality ('_rerun_raw_' = F) uses thermodynamic models and fitness scores from DMS experiments (already processed with MoCHI and DiMSum respectively; see [Required Data](#required-data)) to reproduce all figures in the publication.

# Pipeline Stages

The top-level function **doubledeepms()** is the recommended entry point to the pipeline and reproduces the figures and results from the computational analyses described in the following publication: "A comprehensive map of allosteric sites in a protein domain". See see [Required Data](#required-data) for instructions on how to obtain all required data and miscellaneous files before running the pipeline.

## Stage 1: Evaluate thermodynamic model results

This stage ('doubledeepms_thermo_model_results') evaluates thermodynamic model results and performance comparing to literature free energies.



