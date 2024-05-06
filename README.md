# TotalSeq Analysis Methods For Healthy Control Vs BCR Uveitis Patients
General scripts/code/pipeline/documentation used for analyzing the BCR uveitis single-cell RNA-Seq data

- Authors: Vijay Nagarajan PhD, Pulak Nath PhD
- Affiliation: Laboratory of Immunology, NEI/NIH
- Contact: nagarajanv@nih.gov , pulak.nath@nih.gov
- Description: This workflow takes the raw fastq files, generates the count matrix and performs relevant downstream analysis
- Platform: The workflow BASH commands and R scripts were developed to run in the NIH Biowulf cluster computing facility, but could be reused/reproduced in any linux environment with appropriate changes
-------------------------

### 1. Project Folder Structure

- index.html : (https://pulaknath.github.io/bcr-uveitis/index.html) describing the workflow, with detailed commands, system output and results for the two sample demo workflow analysis
- [Protocols](Protocols) : Contains experiment/protocol details
- [docs](docs) : Contains all the html files for the workflow documentation
- [Results](Results) : Contains all the results for the example data run
- [Tools](Tools) : Contains all the scripts and reference files
- [SampleInformation.tab](SampleInformation.tab) : Contains the sample information for the example data


  [407cfcf7]: https://effective-robot-f7fafb1a.pages.github.io/ "Uveitis TotalSeq Workflow"


### 2. Workflow
The workflow consists of the following components;

1. [Primary analysis](https://pulaknath.github.io/bcr-uveitis/primary-analysis.html)
2. [Quality Control](https://pulaknath.github.io/bcr-uveitis/qualityControl.html)
3. [Normalization](https://pulaknath.github.io/bcr-uveitis/normalization.html)
4. [Clustering](https://pulaknath.github.io/bcr-uveitis/clustering.html)
5. [Differential Expression Analysis](https://pulaknath.github.io/bcr-uveitis/differential.html)
6. [Annotation](https://pulaknath.github.io/bcr-uveitis/annotation.html)
7. [Subcluster analysis](https://github.com/PulakNath/bcr-uveitis/blob/main/Tools/subcluster.R)

### 3. System requirements
- The workflow was developed, tested and run in CentOS Linux 7 (Core) machine, with R version 4.2.2 (2022-10-31), Python 3.8 and Pandas 1.2.4.
- All bash scripts could be run on any unix like operating systems (mac or linux).
- All R scripts could be run on any operating system with R and/or Rstudio and necessary packages installed (all dependencies are provided in the documentation)
- The dependencies for the other included open source applications (SCSA and SCTYPE) are provided in the respective install folders.
- The dependency information for our workflow R scripts are generated using the sessionInfo() command and provided, with specific version numbers, with every corresponding html output.
- Any dependency for our workflow bash scripts are documented in the respective scripts.
- Batch jobs were run using the NIH Biowulf Swarm scripts, but commands for individual jobs are provided in the swarm text file for use with serial job submissions.
- No non-standard hardware is required to run the workflow components.

### 4. Installation guide
The workflow components could be dowloaded and installed from the below link in less than an hour or so depending on the download speed of the user;
https://github.com/PulakNath/bcr-uveitis/archive/refs/heads/main.zip
The workflow components take several hours to run, depending on the number of samples and number of cells in the dataset. Run time for individual components, for the demo data analysis is provided as part of the result/documentation. Detailed setup and any required install guidance is also provided in the workflow documentation.

### 5. Demo data and instructions for use
- All of the raw data, including the two example samples used in the demo documentation are available at the following SRA link;
  - https://www.ncbi.nlm.nih.gov/bioproject/PRJNA855114
- The cellranger processed data (count matrix and metric summary), are available for the two demo example dataset used in the documentation, in the [Results/cellranger](Results/cellranger) folder.
- Detailed instructions to use the workflow for analyzing the demo data is provided in the included documentation, along with the expected output and expected run times.
- 
### 6. How to cite this material
10.5281/zenodo.11117731

### 7. Software license
The license is provided under the [NEI Software Distribution Agreement](LICENSE).

This work utilized the computational resources of the NIH HPC Biowulf cluster.
