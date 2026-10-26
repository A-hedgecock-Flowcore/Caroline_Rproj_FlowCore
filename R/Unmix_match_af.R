#' Unmix a patient sample using matched autofluorescence spectra
#'
#' Extracts the patient ID from the sample filename, finds the corresponding
#' autofluorescence spectra, reports the files being used, and runs AutoSpectral
#' unmixing.
#'
#' @param sample Character string containing the path to one patient `.fcs` file.
#'
#' @return Unmixed .fcs files via  `AutoSpectral::unmix.fcs()`.
#'
#' @details
#' The sample filename must contain a patient identifier in the format
#' `PT` followed by one or more digits, such as `PT160`.
#'
#' `af_matches` must be a data frame or tibble containing the columns
#' `patient_id` and `af_file`.
#'
#' #' The function also requires the objects `spectra`, `asp`, `flow_control`,
#' `fluorophore_variants`, and `unmixed_folder` to be available in the
#' calling environment.
#'
unmix_match_af <- function(sample) {
  patient_id_sample <- str_extract(
    basename(sample),
    "PT[0-9]+"
  )

  if (is.na(patient_id_sample)) {
    stop(
      "Could not identify patient ID in sample: ",
      basename(sample)
    )
  }

  af_match <- af_matches |>
    dplyr::filter(patient_id == patient_id_sample)

  if (nrow(af_match) == 0) {
    stop("No AF spectra found for patient: ", patient_id_sample)
  }

  if (nrow(af_match) > 1) {
    stop("Multiple AF spectra found for patient: ", patient_id_sample)
  }

  af_file_used <- af_match$af_file[[1]]

  cli::cli_alert_success(
    paste("Unmixing sample:", basename(sample))
  )

  cli::cli_alert_success(
    paste("Using AF file:", basename(af_file_used))
  )

  af.spectra <- read.spectra(
    spectra.file = basename(af_file_used),
    spectra.dir = table_spectra_folder
  )

  AutoSpectral::unmix.fcs(
    fcs.file = sample,
    spectra = spectra,
    asp = asp,
    flow.control = flow_control,
    method = "AutoSpectral",
    af.spectra = af.spectra,
    spectra.variants = fluorophore_variants,
    speed = "slow",
    output.dir = unmixed_folder,
    verbose = FALSE,
    pipeline = "joint",
    estimate.time = FALSE
  )
}
