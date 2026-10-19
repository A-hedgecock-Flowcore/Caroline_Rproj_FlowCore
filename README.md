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

