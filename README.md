# Caroline_Rproj_FlowCore

This is a project to create an R script for Caroline Stockwell for a semi-automated analysis pipeline. 
Stated goals for this script are:

## Goals
1. Unmix raw .fcs files via AutoSpectral 
2. Clean unmixed files with semi-automated gating
3. Perform clustering and dimenstion reduction algorithms to aid manual analysis of the the unmixed files
4. Produce an output of:
    1. unimixed files
    2. report of unmixing results
    3. report of clean up gating results
    4. report of clustering and dimReduction steps
    5. final report 


## Set up
The first time this script is run, install the envrionment onto computer:

```{r}
install.packages("renv")
renv::restore()
```

Before running the pipeline, create this file structure:

Input/
├── SCC/
└── Samples/
└── app.r # this is the Shiny app that helps with the AutoSpectral control file

AnalysisInfo/
├── parameters.yml
└── gatingTemplate.csv

R/
├── QCUnmix.R
└── Unmix.R
└── AnalysisPipeline.R

Output/
├── Unmix_files

If the folder names are different than above, update the `parameters.yml` file. 

The pipeline will create the required folders under Output/.