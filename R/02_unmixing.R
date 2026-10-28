# AutoSpectral Variables
asp <- get.autospectral.param(cytometer = "aurora")
control.dir <- single_color_control_folder
control_filename <- list.files(
  path = analysisInfo_folder,
  pattern = analysis_parameters$AutoSpectral_control_file_name,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(control_filename) != 1) {
  stop(
    "Expected exactly one AutoSpectral control file matching pattern '",
    analysis_parameters$AutoSpectral_control_file_name,
    "'. Found ",
    length(control_filename),
    " file(s): ",
    paste(basename(control_filename), collapse = ", ")
  )
}

table_spectra_folder <- file.path(output_folder, "table_spectra")
existing_spectra_files <- list.files(
  table_spectra_folder,
  pattern = "^Automated_autospectral_spectra_[0-9]{8}_[0-9]{6}\\.csv$",
  full.names = TRUE
)
variants_folder <- file.path(output_folder, "variants")

existing_variant_RDS <- list.files(
  path = variants_folder,
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

# Validate patient samples have matching unstained sample
patient_ids_from_samples <- stringr::str_extract(
  basename(Patient_samples),
  analysis_parameters$patient_id_pattern
) |>
  stats::na.omit() |>
  unique()

patient_ids_from_af <- stringr::str_extract(
  basename(Patient_unstained_samples),
  analysis_parameters$patient_id_pattern
) |>
  stats::na.omit() |>
  unique()

missing_af_patients <- setdiff(
  patient_ids_from_samples,
  patient_ids_from_af
)

if (length(missing_af_patients) > 0) {
  warning(
    "Some patient samples do not have matching unstained samples. ",
    "Those samples will be flagged and skipped: ",
    paste(missing_af_patients, collapse = ", ")
  )
}

message(
  "Validated patient autofluorescence samples for: ",
  paste(patient_ids_from_samples, collapse = ", ")
)

# Extract AF signatures from patient samples
if (length(Patient_unstained_samples) == 0) {
  af_matches <- tibble::tibble(
    patient_id = character(),
    af_file = character(),
    af_status = character(),
    af_reason = character()
  )
} else {
  af_matches <- purrr::map_dfr(
    Patient_unstained_samples,
    extract_af_if_missing
  )
}

af_status_file <- file.path(
  output_folder,
  "patient_af_status.csv"
)

readr::write_csv(
  af_matches,
  af_status_file
)

message(
  "Patient AF summary: ",
  sum(af_matches$af_status == "available"),
  " available; ",
  sum(af_matches$af_status == "flagged"),
  " flagged."
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
    variants_folder,
    get.spectral.variants(
      control.dir = normalizePath(
        control.dir,
        winslash = "/",
        mustWork = TRUE
      ),
      control.def.file = normalizePath(
        control_filename,
        winslash = "/",
        mustWork = TRUE
      ),
      asp = asp,
      spectra = spectra,
      output.dir = ".",
      verbose = FALSE
    )
  )
  cli::cli_alert_success("Created new fluorophore variants file.")
}

# Unmix all patient samples with match AF and fluorophore vairants
patient_unmix_status <- purrr::map_dfr(
  Patient_samples,
  unmix_match_af
)

## save unmixing status results for QC report
patient_unmix_status_file <- file.path(
  output_folder,
  "patient_unmixing_status.csv"
)
readr::write_csv(
  patient_unmix_status,
  patient_unmix_status_file
)

## unmixing patient run summary
message(
  "Patient unmixing summary: ",
  sum(patient_unmix_status$status == "completed"),
  " completed; ",
  sum(patient_unmix_status$status == "flagged"),
  " flagged; ",
  sum(patient_unmix_status$status == "failed"),
  " failed."
)
if (!any(patient_unmix_status$status == "completed")) {
  warning(
    "No patient samples were successfully unmixed. ",
    "Review: ",
    patient_unmix_status_file
  )
}

# Unmix any additional samples that are not patient samples
unmixed_files <- list.files(
  path = unmixed_folder,
  full.names = FALSE
)

unmixed_files_names <- paste0(
  tools::file_path_sans_ext(basename(nonpatient_samples)),
  " FastPoisson.fcs" # Change if method below changes
)

samples_to_unmix <- nonpatient_samples[
  !unmixed_files_names %in% unmixed_files
]

if (length(samples_to_unmix) == 0) {
  nonpatient_unmix_status <- tibble::tibble(
    sample = character(),
    status = character(),
    reason = character()
  )
} else {
  nonpatient_unmix_status <- purrr::map_dfr(
    samples_to_unmix,
    function(sample) {
      sample_name <- basename(sample)

      unmix_error <- tryCatch(
        {
          AutoSpectral::unmix.fcs(
            fcs.file = sample,
            spectra = spectra,
            asp = asp,
            flow.control = flow_control,
            method = "FastPoisson",
            output.dir = unmixed_folder,
            verbose = FALSE,
            estimate.time = FALSE
          )

          NULL
        },
        error = function(e) {
          e
        }
      )

      if (is.null(unmix_error)) {
        tibble::tibble(
          sample = sample_name,
          status = "completed",
          reason = NA_character_
        )
      } else {
        warning(
          "FastPoisson unmixing failed for ",
          sample_name,
          ": ",
          conditionMessage(unmix_error)
        )

        tibble::tibble(
          sample = sample_name,
          status = "failed",
          reason = conditionMessage(unmix_error)
        )
      }
    }
  )
}

nonpatient_unmix_status_file <- file.path(
  output_folder,
  "nonpatient_unmixing_status.csv"
)

readr::write_csv(
  nonpatient_unmix_status,
  nonpatient_unmix_status_file
)
