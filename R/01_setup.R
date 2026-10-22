# Check required packages are available
required_packages <- c(
  "tidyverse",
  "yaml",
  "AutoSpectral",
  "AutoSpectralRcpp",
  "Luciernaga",
  "openCyto",
  "flowWorkspace",
  "flowGate",
  "CompensAID"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages. Run renv::restore(). Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

# load necessary packages
library(tidyverse)
library(AutoSpectral)
library(AutoSpectralRcpp)
library(Luciernaga)
library(openCyto)
library(flowWorkspace)
library(flowGate)
library(CompensAID)

# Pull R functions

# file paths
analysisInfo_folder <- file.path("AnalysisInfo")
parameters_file <- list.files(
  analysisInfo_folder,
  "\\.yml$",
  full.names = TRUE
)
analysis_parameters <- yaml::read_yaml(parameters_file)

input_folder <- file.path(analysis_parameters$input_folder)
output_folder <- file.path(analysis_parameters$output_folder)
single_color_control_folder <- file.path(
  input_folder,
  analysis_parameters$single_color_controls_folder
)
samples_folder <- file.path(input_folder, analysis_parameters$samples_folder)
unmixed_folder <- file.path(output_folder, analysis_parameters$unmixed_folder)

Patient_unstained_samples <- list.files(
  single_color_control_folder,
  pattern = analysis_parameters$patient_sample_pattern,
  full.names = TRUE,
  ignore.case = TRUE
)

Patient_samples <- list.files(
  path = samples_folder,
  pattern = analysis_parameters$patient_sample_pattern,
  full.names = TRUE,
  ignore.case = TRUE
)

CleanUpGatingTemplateFile <- list.files(
  path = analysisInfo_folder,
  pattern = analysis_parameters$gatingTemplate,
  full.names = TRUE,
  ignore.case = TRUE
)
