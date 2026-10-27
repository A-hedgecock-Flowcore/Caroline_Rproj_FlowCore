#' Unmix one patient sample using matched autofluorescence spectra
#'
#' Extracts the patient identifier from the sample filename, finds the
#' corresponding autofluorescence spectra in `af_matches`, and runs
#' `AutoSpectral::unmix.fcs()` when exactly one matching AF file is available.
#'
#' Missing patient IDs, missing AF spectra, duplicate AF matches, AF-reading
#' errors, and AutoSpectral errors are returned as status records rather than
#' stopping the entire pipeline. This allows other patient samples to
#' continue processing.
#'
#' The function depends on the following objects being available in its
#' calling environment:
#' \describe{
#'   \item{af_matches}{A data frame or tibble containing `patient_id` and
#'   `af_file` columns.}
#'   \item{table_spectra_folder}{Directory containing AF spectra files.}
#'   \item{spectra}{The AutoSpectral reference spectra.}
#'   \item{asp}{An AutoSpectral parameter object.}
#'   \item{flow_control}{The AutoSpectral flow-control object.}
#'   \item{fluorophore_variants}{The AutoSpectral fluorophore variant object.}
#'   \item{unmixed_folder}{Directory where unmixed FCS files are written.}
#' }
#'
#' Patient identifiers must occur in the sample filename in the form
#' `"PT"` followed by one or more digits, such as `"PT169"`.
#'
#' @param sample Character string containing the path to one patient FCS
#'   sample.
#'
#' @return A one-row tibble with the following columns:
#' \describe{
#'   \item{sample}{The basename of the patient sample.}
#'   \item{patient_id}{The extracted patient identifier.}
#'   \item{status}{One of `"completed"`, `"flagged"`, or `"failed"`.}
#'   \item{reason}{A description of the issue, or `NA` when unmixing
#'   completed successfully.}
#'   \item{af_file}{The AF spectra filename used, or `NA` if no file was
#'   available.}
#' }
#'
#' @details
#' A sample is processed as follows:
#' \enumerate{
#'   \item The patient identifier is extracted from the sample filename.
#'   \item Matching AF records are found in `af_matches`.
#'   \item Exactly one AF spectra file is read.
#'   \item `AutoSpectral::unmix.fcs()` is run using the matched AF spectra.
#'   \item A status tibble is returned for use in QC summaries.
#' }
#'
#' Samples with no matching AF spectra or multiple matching AF records are
#' flagged and skipped. Samples for which AF reading or AutoSpectral unmixing
#' fails are returned with status `"failed"`.
#'
#' @examples
#' \dontrun{
#' result <- unmix_match_af(
#'   sample = "Input/Samples/11 PT169 immune panel.fcs"
#' )
#'
#' result
#' }
#'
#' @export
unmix_match_af <- function(sample) {
  sample_name <- basename(sample)

  patient_id_sample <- stringr::str_extract(
    sample_name,
    analysis_parameters$patient_id_pattern
  )

  make_status <- function(
    status,
    reason = NA_character_,
    af_file = NA_character_
  ) {
    tibble::tibble(
      sample = sample_name,
      patient_id = patient_id_sample,
      status = status,
      reason = reason,
      af_file = af_file
    )
  }

  if (is.na(patient_id_sample)) {
    warning(
      "Could not identify patient ID in sample: ",
      sample_name
    )

    return(
      make_status(
        status = "flagged",
        reason = "Could not identify patient ID"
      )
    )
  }

  af_match <- af_matches |>
    dplyr::filter(
      patient_id == patient_id_sample
    )

  if (nrow(af_match) == 0) {
    reason <- "No matching patient AF record"

    warning(
      reason,
      " for ",
      patient_id_sample,
      ". Skipping sample: ",
      sample_name
    )

    return(
      make_status(
        status = "flagged",
        reason = reason
      )
    )
  }

  if (nrow(af_match) > 1) {
    reason <- "Multiple matching patient AF records"

    warning(
      reason,
      " for ",
      patient_id_sample,
      ". Skipping sample: ",
      sample_name
    )

    return(
      make_status(
        status = "flagged",
        reason = reason
      )
    )
  }

  if (
    is.na(af_match$af_file[[1]]) ||
      !identical(af_match$af_status[[1]], "available")
  ) {
    reason <- af_match$af_reason[[1]]

    if (is.na(reason) || !nzchar(reason)) {
      reason <- "Patient AF spectra are unavailable"
    }

    warning(
      reason,
      " Skipping sample: ",
      sample_name
    )

    return(
      make_status(
        status = "flagged",
        reason = reason
      )
    )
  }

  af_file_used <- af_match$af_file[[1]]

  message(
    "Unmixing patient sample: ",
    sample_name
  )

  message(
    "Using AF file: ",
    basename(af_file_used)
  )

  af_spectra <- tryCatch(
    {
      AutoSpectral::read.spectra(
        spectra.file = basename(af_file_used),
        spectra.dir = table_spectra_folder
      )
    },
    error = function(e) {
      e
    }
  )

  if (inherits(af_spectra, "error")) {
    reason <- paste(
      "Could not read AF spectra:",
      conditionMessage(af_spectra)
    )

    warning(
      reason,
      " Sample: ",
      sample_name
    )

    return(
      make_status(
        status = "failed",
        reason = reason,
        af_file = basename(af_file_used)
      )
    )
  }

  unmix_error <- tryCatch(
    {
      AutoSpectral::unmix.fcs(
        fcs.file = sample,
        spectra = spectra,
        asp = asp,
        flow.control = flow_control,
        method = "AutoSpectral",
        af.spectra = af_spectra,
        spectra.variants = fluorophore_variants,
        speed = "slow",
        output.dir = unmixed_folder,
        verbose = FALSE,
        pipeline = "joint",
        estimate.time = FALSE
      )

      NULL
    },
    error = function(e) {
      e
    }
  )

  if (!is.null(unmix_error)) {
    reason <- paste(
      "AutoSpectral unmixing failed:",
      conditionMessage(unmix_error)
    )

    warning(
      reason,
      " Sample: ",
      sample_name
    )

    return(
      make_status(
        status = "failed",
        reason = reason,
        af_file = basename(af_file_used)
      )
    )
  }

  make_status(
    status = "completed",
    reason = NA_character_,
    af_file = basename(af_file_used)
  )
}
