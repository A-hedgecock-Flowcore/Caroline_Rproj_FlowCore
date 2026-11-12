# Clustering and Dimenstion Reduction steps after QC of data has been verified

# Setup
library(flowWorkspace)
library(Spectre)
library(tidyverse)
library(data.table)

## file.paths
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

output_folder <- file.path("Output")
analysis_outputs <- file.path(output_folder, "analysis_outputs")
if (!dir.exists(analysis_outputs)) {
  dir.create(analysis_outputs)
}

# choose latest gatingSet and load
output_folder_contents <- list.dirs(
  output_folder,
  full.names = TRUE,
  recursive = FALSE
)
PT_gatingSets <- output_folder_contents[str_detect(
  output_folder_contents,
  pattern = "PT_GatingSet"
)]

if (length(PT_gatingSets) > 1) {
  PT_gs_path <- sort(PT_gatingSets)[length(PT_gatingSets)]
  cli::cli_alert_success(paste0("Using GatingSet ", basename(PT_gs_path)))
} else {
  PT_gs_path <- PT_gatingSets
}

PT_gs <- load_gs(PT_gs_path)

# Transfer gs to a data.table
cleaned_data <- gs_pop_get_data(PT_gs, "singlets", inverse.transform = FALSE)
samples_names <- sampleNames(PT_gs)

cleaned_data_list <- map(samples_names, \(sample_id) {
  pop_data <- exprs(cleaned_data[[sample_id]])
  dt <- as.data.table(pop_data)
  dt[, Sample := sample_id]
  return(dt)
})

cleaned_data <- rbindlist(cleaned_data_list, use.names = TRUE, fill = TRUE)
cleaned_data[, Sample := str_remove(Sample, " AutoSpectral\\.fcs$")]

message(
  "Cleaned data from samples ",
  paste(unique(cleaned_data$Sample), collapse = ", "),
  "."
)
