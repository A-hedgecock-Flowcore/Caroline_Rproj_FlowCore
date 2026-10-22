# AutoSpectral Variables
asp <- get.autospectral.param(cytometer = "aurora")
control.dir <- single_color_control_folder
control_filename <- list.files(
  path = output_folder,
  pattern = analysis_parameters$AutoSpectral_control_file_name,
  full.names = TRUE,
  ignore.case = TRUE
)
table_spectra_folder <- file.path(output_folder, "table_spectra")
existing_spectra_files <- list.files(
  table_spectra_folder,
  pattern = "^Automated_autospectral_spectra_[0-9]{8}_[0-9]{6}\\.csv$",
  full.names = TRUE
)

# R helper functions
source(file.path("R", "latest_af_file.R"))
source(file.path("R", "spectra_extraction.R"))
source(file.path("R", "extract_af_if_missing.R"))
source(file.path("R", "unmix_match_af.R"))


# Validate control file
issues <- check.control.file(
  control.dir = control.dir,
  control.def.file = control_filename,
  asp,
  verbose = FALSE
)

# Spectral Extraction

spectra <- spectra_extraction(
  output_folder = output_folder,
  control_dir = control.dir,
  control_filename = control_filename,
  asp = asp
)

# load in the AutoSpectral Control file made during set up (see README)
flow_control <- reload.flow.control(control.dir, control_filename, asp)

# Extract AF signatures from patient samples
af_matches <- map_dfr(
  Patient_unstained_samples,
  extract_af_if_missing
)
