#' Unmix a patient sample using matched autofluorescence spectra
#'
#' Extracts the patient ID from the sample filename, finds the corresponding
#' autofluorescence spectra, reports the files being used, and runs AutoSpectral
#' unmixing.
#'
#' @param sample Character string containing the path to one patient `.fcs` file.
#'
#' @return The result returned by `AutoSpectral::unmix.fcs()`.
#'
#' @details
#' The sample filename must contain a patient identifier in the format
#' `PT` followed by one or more digits, such as `PT160`.
#'
#' The object `patient_af` must be a named list of autofluorescence spectra,
#' and `af_files` must be a named vector of AF filenames. Both must use
#' patient IDs as names.
#'
#' @examples
#' \dontrun{
#' unmixTest(Patient_samples[3])
#' }
#'
Unmix_match_af <- function(sample) {
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

  if (!patient_id_sample %in% names(patient_af)) {
    stop("No AF spectra found for patient: ", patient_id_sample)
  }

  af_for_patient <- patient_af[[patient_id_sample]]
  af_file_used <- af_files[[patient_id_sample]]

  message("Unmixing sample: ", basename(sample))
  message("Using AF file: ", basename(af_file_used))

  cli::cli_alert_success(
    paste("Unmixing sample:", basename(sample))
  )

  cli::cli_alert_success(
    paste("Using AF file:", basename(af_file_used))
  )

  AutoSpectral::unmix.fcs(
    fcs.file = sample,
    spectra = spectra,
    asp = asp,
    flow.control = flow_control,
    method = "AutoSpectral",
    af.spectra = af_for_patient,
    spectra.variants = fluorophore_variants,
    speed = "slow",
    output.dir = unmixed_folder,
    verbose = FALSE,
    pipeline = "joint",
    estimate.time = FALSE
  )
}
