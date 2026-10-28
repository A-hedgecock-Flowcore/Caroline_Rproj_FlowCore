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

# file paths
analysisInfo_folder <- file.path("AnalysisInfo")
parameters_file <- list.files(
  analysisInfo_folder,
  "\\.yml$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(parameters_file) != 1) {
  stop(
    "Expected exactly one (1) .yml (YAML) parameter file in ",
    analysisInfo_folder,
    ". Found ",
    length(parameters_file),
    "."
  )
}

analysis_parameters <- yaml::read_yaml(parameters_file)

input_folder <- file.path(analysis_parameters$input_folder)

output_folder <- file.path(analysis_parameters$output_folder)

single_color_control_folder <- file.path(
  input_folder,
  analysis_parameters$single_color_controls_folder
)

samples_folder <- file.path(
  input_folder,
  analysis_parameters$samples_folder
)

unmixed_folder <- file.path(
  output_folder,
  analysis_parameters$unmixed_folder
)

required_input_dirs <- c(
  single_color_control_folder,
  samples_folder
)

missing_input_dirs <- required_input_dirs[!dir.exists(required_input_dirs)]

if (length(missing_input_dirs) > 0) {
  stop(
    "Missing input directories: ",
    paste(missing_input_dirs, collapse = ", ")
  )
}

required_output_directories <- c(
  output_folder,
  unmixed_folder,
  file.path(output_folder, "table_spectra"),
  file.path(output_folder, "compensAID"),
  file.path(output_folder, "logs"),
  file.path(output_folder, "variants")
)

purrr::walk(
  required_output_directories,
  dir.create,
  showWarnings = FALSE,
  recursive = TRUE
)

Patient_unstained_samples <- list.files(
  single_color_control_folder,
  pattern = analysis_parameters$patient_unstained_pattern,
  full.names = TRUE,
  ignore.case = TRUE
)

Patient_samples <- list.files(
  path = samples_folder,
  pattern = analysis_parameters$patient_sample_pattern,
  full.names = TRUE,
  ignore.case = TRUE
)

All_samples <- list.files(
  samples_folder,
  pattern = "\\.fcs$",
  full.names = TRUE,
  ignore.case = TRUE
)

nonpatient_samples <- All_samples[
  !str_detect(
    basename(All_samples),
    analysis_parameters$patient_id_pattern
  )
]

CleanUpGatingTemplateFile <- list.files(
  path = analysisInfo_folder,
  pattern = analysis_parameters$gatingTemplate,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(CleanUpGatingTemplateFile) != 1) {
  stop(
    "Expected exactly one gating template file. Found ",
    length(CleanUpGatingTemplateFile),
    "."
  )
}

if (length(Patient_samples) == 0) {
  warning(
    "No patient samples were found in: ",
    samples_folder
  )
}

if (length(All_samples) == 0) {
  stop(
    "No .fcs files were found in: ",
    samples_folder
  )
}

message(
  "Setup complete: ",
  length(Patient_samples),
  " patient samples, ",
  length(Patient_unstained_samples),
  " patient unstained samples, and ",
  length(nonpatient_samples),
  " nonpatient samples."
)
