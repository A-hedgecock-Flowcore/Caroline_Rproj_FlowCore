#' Extract or identify a patient's autofluorescence spectra
#'
#' Checks whether an autofluorescence spectra file already exists for the
#' patient represented by an unstained sample. If a matching file exists,
#' the most recently modified file is used. If no matching file exists,
#' autofluorescence spectra are extracted using `AutoSpectral::get.af.spectra()`.
#'
#' @param sample Character string containing the path to one unstained
#'   patient sample.
#'
#' @return A tibble with one row and two columns:
#' \describe{
#'   \item{patient_id}{The patient identifier extracted from the sample
#'   filename.}
#'   \item{af_file}{The path to the AF spectra file used or created.}
#' }
#'
#' @details
#' The sample filename must contain a patient identifier in the format
#' `"PT"` followed by one or more digits, such as `"PT160"`.
#'
#' This function depends on the following objects being available in the
#' calling environment:
#' \describe{
#'   \item{table_spectra_folder}{Directory for AF spectra CSV files.}
#'   \item{asp}{AutoSpectral object.}
#'   \item{spectra}{AutoSpectral spectra matrix.}
#'   \item{plot.dir}{Directory for generated AF plots.}
#' }
#'

extract_af_if_missing <- function(sample) {
  patient_id <- str_extract(
    basename(sample),
    analysis_parameters$patient_id_pattern
  )

  if (is.na(patient_id)) {
    stop("Could not identify patient ID in : ", basename(sample))
  }

  af_file <- latest_af_file(
    patient_id = patient_id,
    table_spectra_folder = table_spectra_folder
  )

  if (length(af_file) > 0) {
    message(
      "Using exisiting AF file for ",
      patient_id,
      ": ",
      basename(af_file)
    )
  } else {
    message(
      "No AF file found for ",
      patient_id,
      ". Extracting AF spectra."
    )

    plot.dir <- file.path(output_folder, "figure_autofluorscence")

    get.af.spectra(
      unstained.sample = sample,
      asp = asp,
      spectra = spectra,
      plot.dir = plot.dir,
      table.dir = table_spectra_folder,
      verbose = FALSE
    )

    # find newly created af file
    af_file <- latest_af_file(
      patient_id = patient_id,
      table_spectra_folder = table_spectra_folder
    )

    if (length(af_file) == 0) {
      stop(
        "AF extraction completed, but no AF file found for :",
        patient_id
      )
    }
    message(
      "Created AF file for ",
      patient_id,
      ": ",
      basename(af_file)
    )
  }

  tibble::tibble(
    patient_id = patient_id,
    af_file = af_file
  )
}
