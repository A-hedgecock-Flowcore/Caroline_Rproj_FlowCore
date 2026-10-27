#' Extract or identify a patient's autofluorescence spectra
#'
#' Finds an existing autofluorescence spectra file for a patient unstained
#' sample or extracts one using `AutoSpectral::get.af.spectra()`.
#'
#' AF extraction problems are returned as status records rather than stopping
#' the entire pipeline. This allows other patient samples to continue.
#'
#' @param sample Character string containing the path to one patient
#'   unstained FCS sample.
#'
#' @return A one-row tibble with:
#' \describe{
#'   \item{patient_id}{The extracted patient identifier.}
#'   \item{af_file}{The AF spectra file path, or `NA` if unavailable.}
#'   \item{af_status}{One of `"available"` or `"flagged"`.}
#'   \item{af_reason}{Description of the issue, or `NA` if successful.}
#' }
#'
#' @details
#' The sample filename must contain a patient identifier matching
#' `analysis_parameters$patient_id_pattern`, such as `"PT169"`.
#'
#' The function depends on the following objects being available in the
#' calling environment:
#' \describe{
#'   \item{analysis_parameters}{The YAML configuration list.}
#'   \item{table_spectra_folder}{Directory containing AF spectra files.}
#'   \item{output_folder}{Main output directory.}
#'   \item{asp}{The AutoSpectral parameter object.}
#'   \item{spectra}{The AutoSpectral reference spectra.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- extract_af_if_missing(
#'   sample = "Input/SCC/PT169 unstained (Cells).fcs"
#' )
#'
#' result
#' }
#'
#' @export
extract_af_if_missing <- function(sample) {
  sample_name <- basename(sample)

  patient_id <- stringr::str_extract(
    sample_name,
    analysis_parameters$patient_id_pattern
  )

  make_status <- function(
    patient_id,
    af_file = NA_character_,
    af_status = "flagged",
    af_reason = NA_character_
  ) {
    tibble::tibble(
      patient_id = patient_id,
      af_file = af_file,
      af_status = af_status,
      af_reason = af_reason
    )
  }

  if (is.na(patient_id)) {
    reason <- paste(
      "Could not identify patient ID in unstained sample:",
      sample_name
    )

    warning(reason)

    return(
      make_status(
        patient_id = NA_character_,
        af_reason = reason
      )
    )
  }

  af_file <- latest_af_file(
    patient_id = patient_id,
    table_spectra_folder = table_spectra_folder
  )

  if (length(af_file) > 0) {
    message(
      "Using existing AF file for ",
      patient_id,
      ": ",
      basename(af_file)
    )

    return(
      make_status(
        patient_id = patient_id,
        af_file = af_file,
        af_status = "available"
      )
    )
  }

  message(
    "No AF file found for ",
    patient_id,
    ". Extracting AF spectra."
  )

  plot_dir <- file.path(
    output_folder,
    "figure_autofluorescence"
  )

  dir.create(
    plot_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  extraction_error <- tryCatch(
    {
      AutoSpectral::get.af.spectra(
        unstained.sample = sample,
        asp = asp,
        spectra = spectra,
        plot.dir = plot_dir,
        table.dir = table_spectra_folder,
        verbose = FALSE
      )

      NULL
    },
    error = function(e) {
      e
    }
  )

  if (!is.null(extraction_error)) {
    reason <- paste(
      "AF extraction failed:",
      conditionMessage(extraction_error)
    )

    warning(
      reason,
      " Sample: ",
      sample_name
    )

    return(
      make_status(
        patient_id = patient_id,
        af_reason = reason
      )
    )
  }

  af_file <- latest_af_file(
    patient_id = patient_id,
    table_spectra_folder = table_spectra_folder
  )

  if (length(af_file) == 0) {
    reason <- paste(
      "AF extraction completed, but no AF spectra file was found for:",
      patient_id
    )

    warning(reason)

    return(
      make_status(
        patient_id = patient_id,
        af_reason = reason
      )
    )
  }

  message(
    "Created AF file for ",
    patient_id,
    ": ",
    basename(af_file)
  )

  make_status(
    patient_id = patient_id,
    af_file = af_file,
    af_status = "available"
  )
}
