# Spectral Cytometry settings
asp <- get.autospectral.param(cytometer = "aurora")
control.dir <- single_color_control_folder
control_filename <- list.files(
  path = output_folder,
  pattern = analysis_parameters$AutoSpectral_control_file_name,
  full.names = TRUE,
  ignore.case = TRUE
)

# Validate control file
issues <- check.control.file(
  control.dir = control.dir,
  control.def.file = control_filename,
  asp
)

# Spectral Extraction
table_spectra_folder <- file.path(output_folder, "table_spectra")
existing_spectra_files <- list.files(
  table_spectra_folder,
  pattern = "^Automated_autospectral_spectra_[0-9]{8}_[0-9]{6}\\.csv$",
  full.names = TRUE
)

# does a spectra file already exist?
if (base::length(existing_spectra_files) > 0) {
  file_info <- file.info(existing_spectra_files)

  # If several files exist, use the most recently modified one
  latest_file <- rownames(
    file_info[which.max(file_info$mtime), , drop = FALSE]
  )

  spectra <- read.spectra(
    spectra.file = basename(latest_file),
    spectra.dir = dirname(latest_file)
  )
  message("Using spectra file: ", basename(latest_file))
} else {
  # if no file exists, create
  full_control.dir <- normalizePath(control.dir, mustWork = TRUE)
  full_control_filename <- normalizePath(control_filename, mustWork = TRUE)

  spectra <- withr::with_dir(
    output_folder,
    get.spectra.automated(
      control.dir = full_control.dir,
      control.def.file = full_control_filename,
      asp = asp,
      figures = TRUE,
      verbose = FALSE
    )
  )
  message("Created new AutoSpectral spectra.")
}

# load in the AutoSpectral Control file made during set up (see README)
flow_control <- reload.flow.control(control.dir, control_filename, asp)

# Extract AF signatures from patient samples
af_files <- list.files(
  path = table_spectra_folder, # find exisitng af files
  pattern = "Autofluorescence spectra.*\\.csv$",
  full.names = FALSE,
  ignore.case = TRUE
)

plot.dir <- file.path(output_folder, "figure_autofluorscence")

patient_af <- map(
  Patient_unstained_samples,
  \(sample) {
    get.af.spectra(
      unstained.sample = sample,
      asp = asp,
      spectra = spectra,
      plot.dir = plot.dir,
      table.dir = table_spectra_folder,
      verbose = FALSE
    )
  }
)

# Extract fluorophore variants per-cell
existing_variant_RDS <- list.files(
  path = output_folder,
  pattern = "\\.rds$",
  full.names = TRUE
)

if (base::length(existing_variant_RDS) > 0) {
  fluorophore_variants <- read_rds(existing_variant_RDS)
  message("Using fluorophore variant file ", basename(existing_variant_RDS))
} else {
  fluorophore_variants <- get.spectral.variants(
    control.dir = control.dir,
    control.def.file = control_filename,
    asp = asp,
    spectra = spectra,
    output.dir = output_folder
  )
}
