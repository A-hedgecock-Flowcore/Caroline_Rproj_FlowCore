# Caroline_Rproj_FlowCore

This is a project to create an R script for Caroline Stockwell for a semi-automated analysis pipeline. Stated goals for this script are:

## Goals

1.  Unmix raw .fcs files via AutoSpectral
2.  Clean unmixed files with semi-automated gating
3.  Perform clustering and dimension reduction algorithms to aid manual analysis of the the unmixed files
4.  Produce an output of:
    1.  unimixed files
    2.  report of unmixing results
    3.  report of clean up gating results
    4.  report of clustering and dimReduction steps
    5.  final report

## Set up

### Restore environment

The first time this script is run, install the environment onto computer:

```{r}
install.packages("renv")
renv::restore()
```

### File Structure

Before running the pipeline, make sure this file structure is intact:

**FIX WHEN COMPLETE**

**Input**/

├── SCC/ \# single color controls

└── Samples/

└── app.r \# this is the Shiny app that helps with the AutoSpectral control file.

**AnalysisInfo**/

├── parameters.yml

└── gatingTemplate.csv

**R**/

├── 01_setup.R

└── 02_AutoSpectral.R

└── QCUnmix.R

└── Unmix.R

└── AnalysisPipeline.R

**Output**/

├── Unmix_files

└── control file \# see creation below

If the folder names are different than above, update the `parameters.yml` file.

The pipeline will create any additional required folders under Output/.

### AutoSpectral Control File

See Control File information here: [Control File- AutoSpectral](https://drcytometer.github.io/AutoSpectral/articles/02_Control_File_example.html) The main structural elements to the control file are the following columns. Note, we will be using the `get.spectra.automated()` function, so the "legacy" columns are not needed. Read more in linked article.

- **filename**: The actual .fcs filename
- **fluorophore**: The name of the fluorophore (or AF in case of the unstained)
- **marker**: The name of the marker (or left empty in the case of the unstained)
- **channel**: The peak detector channel of the fluorophore (left empty for the unstained).
- **control.type**: Whether the unmixing control is cells or beads (used for matching positive and negatives for background subtraction)
- **universal.negative**: By default empty, if you want to use an universal negative, copy the filename of the corresponding unstained into this spot on a per fluorophore basis
- large.gate: Logical value, by default FALSE, set to TRUE to expand the gating area
- gate.name: Gate name for use in signature retrieval
- gate.define: Logical, permits you to set the general gate region as defined by gate.name
- is.viability: Whether your control is a viability dye.

Edit/use the Shiny app to amend the file as necessary. In the future, we can add programmatic edits if wanted.

**Overall, I suggest the Shiny app when possible as it it much more manageable!**

```{r}
#| label: Shiny Helper App AutoSpectral. Creates and validates the control file. 
shiny::runApp(appDir = input_folder, launch.browser = TRUE)
```

```{r}
#| label: Create Control file manually

control.dir <- file.path("Input", "SCC")
asp <- asp <- AutoSpectral::get.autospectral.param(cytometer = "aurora")
output.dir = file.path("AnalysisInfo")
filename <- "fcs_control_file" # if changed, update parameters.yml file 

create.control.file(
  control.dir = control.dir,
  asp,
  output.dir = output.dir,
  filename = filename
)
```

If fluorophores exist in the `AutoSpectral` database they will be auto-filled in the control file. If not, you will need to manually add it in.

The control file needs be checked for any glaring errors. This will also be run during the `AnalysisPipeline.R` in the `AutoSpectral.R` script.

```{r}
#| label: Manual control file needs to be validated if before moving forward

control.dir <- file.path("Input", "SCC")
asp <- asp <- AutoSpectral::get.autospectral.param(cytometer = "aurora")
filename <- file.path("AnalysisInfo", "fcs_control_file".csv) # if changed, update parameters.yml file 


check.control.file(
  control.dir,
  control.def.file = filename,
  asp
)
```