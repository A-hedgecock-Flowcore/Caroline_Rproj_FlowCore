# Custom Flow Cytometry Analysis Pipeline

This is a project to creates a semi-automated flow cytometry analysis pipeline for Caroline Stockwell. Stated goals for this script are:

## Goals

The pipeline is designed to:

1.  Unmix raw .fcs files via AutoSpectral
2.  Clean unmixed files with semi-automated gating
3.  Perform clustering and dimension reduction algorithms to aid manual analysis of the unmixed files
4.  Produce:
    - unmixed files
    - report of unmixing results
    - report of clean up gating results
    - report of clustering and dimReduction steps
    - final report

## Set up

### Installation

This project uses [`renv`](https://rstudio.github.io/renv/) to manage the R package environment.

The first time the project is run, restore the project environment:

```{r}
install.packages("renv")
renv::restore()
```

Run all commands from the project root directory.

### File Structure

Before running the pipeline, make sure this file structure is intact:

``` text
Project Root/
├── Input/
│   ├── Samples/
│   ├── SCC/
│   └── app.R
├── R/
│   ├── 01_setup.R
│   ├── 02_unmixing.R
│   ├── 03_gating.R
│   ├── 04_compensAID.R
│   ├── 05_export_qc.R
│   ├── extract_af_if_missing.R
│   ├── get_gated_flowFrame.R
│   ├── latest_af_file.R
│   ├── qc_export_helper.R
│   ├── run_compensaid.R
│   ├── spectra_extraction.R
│   └── Unmix_match.R
├── AnalysisInfo/
│   ├── fcs_control_file.csv
│   ├── gatingTemplate.csv
│   └── parameters.yml
├── QCPipeline.R
├── QCDashboard.R
├── QCReport.R
├── index.qmd
├── renv.lock
├── _quarto.yml
└── README.md
```

If the folder names are different than above, update the `parameters.yml` file.

The pipeline will create all additional required folders under `Output/` (created by pipeline).

The R files can all be found in [GitHub](https://github.com/A-hedgecock-Flowcore/Caroline_Rproj_FlowCore).

### AutoSpectral Control File

See Control File information here: [Control File- AutoSpectral](https://drcytometer.github.io/AutoSpectral/articles/02_Control_File_example.html) The main structural elements to the control file are the following columns. Note, we will be using the `get.spectra.automated()` function, so the "legacy" columns are not needed. Read more in linked article.

| Column | Description |
|------------------------|------------------------------------------------|
| `filename` | The actual .fcs filename |
| `fluorophore` | The name of the fluorophore (or AF/negative in case of the unstained) |
| `marker` | The name of the marker (or left empty in the case of the unstained) |
| `channel` | The peak detector channel of the fluorophore (left empty for the unstained). |
| `control.type` | Whether the unmixing control is cells or beads (used for matching positive and negatives for background subtraction) |
| `universal.negative` | By default empty, if you want to use an universal negative, copy the filename of the corresponding unstained into this spot on a per fluorophore basis |
| `large.gate` | Logical value, by default FALSE, set to TRUE to expand the gating area. Deprecated when using `get.spectra.automated()` |
| `gate.name` | Gate name for use in signature retrieval. Deprecated when using `get.spectra.automated()` |
| `gate.define` | Logical, permits you to set the general gate region as defined by gate.name. Deprecated when using `get.spectra.automated()` |
| `is.viability` | Whether your control is a viability dye. Deprecated when using `get.spectra.automated()` |

Edit/use the Shiny app to amend the file as necessary. In the future, we can add programmatic edits if wanted.

If fluorophores exist in the `AutoSpectral` database they will be auto-filled in the control file. If not, you will need to manually add it in.

**Overall, I suggest the Shiny app when possible as it is much more manageable!**

#### Shiny App help for control file

```{r}
#| label: Shiny Helper App AutoSpectral. Creates and validates the control file. 
shiny::runApp(appDir = input_folder, launch.browser = TRUE)
```

#### Manually Create the control file

```{r}
#| label: Create Control file manually

control.dir <- file.path("Input", "SCC")
asp <- AutoSpectral::get.autospectral.param(cytometer = "aurora")
output.dir <- file.path("AnalysisInfo")
filename <- "fcs_control_file" # if changed, update parameters.yml file 

create.control.file(
  control.dir = control.dir,
  asp,
  output.dir = output.dir,
  filename = filename
)
```

The control file needs be checked for any glaring errors. This will also be run during the `QCPipeline.R` in the `02_unmixing.R` script.

```{r}
#| label: Manual control file needs to be validated if before moving forward

control.dir <- file.path("Input", "SCC")
asp <- AutoSpectral::get.autospectral.param(cytometer = "aurora")
filename <- file.path("AnalysisInfo", "fcs_control_file.csv") # if changed, update parameters.yml file 


AutoSpectral::check.control.file(
  control.dir,
  control.def.file = filename,
  asp
)
```

### Parameters.yml

The `AnalysisInfo/parameters.yml` file contains project-specific paths and settings, including the input folders, output folders, and AutoSpectral control-file name. Update this file if the project structure changes.

| Parameter | Description |
|-------------------------|----------------------------------------------|
| `input_folder` | Folder name of the folder the files to run the pipeline on. Default is "Input". |
| `output_folder` | Folder name of the folder to put all outputs of the pipeline into. Default is "Output". |
| `single_color_controls_folder` | Folder name within the `input_folder` that contains the single color controls and unstained controls of the experiment. Default is "SCC" |
| `samples_folder` | Folder name within the `input_folder` that contains all other .fcs files that are **not** single stained controls or unstained controls. Default is "Samples" |
| `unmixed_folder` | Folder name within the `output_folder` that will house the unmixed files after being unmixed via AutoSpectral. Default is "Unmixed_files" |
| `patient_id_pattern` | RegEx pattern to identify patient samples. The default is based on the pattern "PT160" |
| `patient_sample_pattern` | RegEx pattern to find the patient .fcs files. The default is based on the pattern "9 PT160 immune panel.fcs" |
| `patient_unstained_pattern` | RegEx pattern to find the unstained patient .fcs files. The default is based on the pattern "PT160 unstained (Cells).fcs" |
| `gatingTemplate` | Filename of the `openCyto` gatingTemplate .csv file. Default is "gatingTemplate.csv" |
| `AutoSpectral_control_file_name` | Filename of the `AutoSpectral` control file. Default is "fcs_control_file" |
| `AutoSpectral_control_file_validated` | Logical. If the R shiny app was used to create the control file, the file is automatically validated. If not, needs to be validated prior to running unmixing. Default is TRUE. |

## Running the QC Pipeline

Run the pipeline from the project root directory.

```{r}
source("QCPipeline.R")
```
