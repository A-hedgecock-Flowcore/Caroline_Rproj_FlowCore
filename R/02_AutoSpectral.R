# AutoSpectral Variables
asp <- get.autospectral.param(cytometer = "aurora")
control.dir <- single_color_control_folder
control_filename <- list.files(
  path = analysisInfo_folder,
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
existing_variant_RDS <- list.files(
  path = output_folder,
  pattern = "Spectral_variants\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE
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

# Fluorophore Variant extraction
if (base::length(existing_variant_RDS) > 0) {
  latest_variant_file <- existing_variant_RDS[
    which.max(file.info(existing_variant_RDS)$mtime)
  ]

  fluorophore_variants <- read_rds(latest_variant_file)
  message("Using fluorophore variant file ", basename(latest_variant_file))
} else {
  fluorophore_variants <- withr::with_dir(
    output_folder,
    get.spectral.variants(
      control.dir = fs::path_abs(control.dir, start = ".."),
      control.def.file = fs::path_abs(control_filename, start = ".."),
      asp = asp,
      spectra = spectra,
      output.dir = ".",
      verbose = FALSE
    )
  )
  cli::cli_alert_success("Created new fluorophore variants file.")
}

# Unmix all patient samples with match AF and fluorophore vairants
walk(Patient_samples, unmix_match_af)

# Unmix any additional samples that are not patient samples
PT_unmixed_files <- list.files(
  path = unmixed_folder,
  full.names = FALSE
)

samples_to_unmix <- nonpatient_samples[
  !basename(nonpatient_samples) %in% basename(PT_unmixed_files)
]

walk(
  samples_to_unmix,
  \(sample) {
    unmix.fcs(
      fcs.file = sample,
      spectra = spectra,
      asp = asp,
      flow.control = flow_control,
      method = "FastPoisson", # not Autospectral; does not conduct per-cell AF extraction. Double check this is fine for controls
      spectra.variants = fluorophore_variants,
      speed = "slow",
      output.dir = unmixed_folder,
      verbose = FALSE,
      pipeline = "joint",
      estimate.time = FALSE
    )
  }
)
