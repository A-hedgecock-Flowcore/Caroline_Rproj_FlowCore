#' Find the latest autofluorescence spectra file for a patient
#'
#' Searches the spectra output directory for autofluorescence spectra files
#' associated with a specified patient ID. If multiple matching files exist,
#' the most recently modified file is returned.
#'
#' @param patient_id Character string containing the patient identifier,
#'   such as `"PT160"`.
#' @param table_spectra_folder Character string specifying the directory
#'   containing the autofluorescence spectra CSV files.
#'
#' @return A character string containing the path to the most recently
#'   modified matching AF file. Returns `character(0)` if no matching file
#'   is found.
#'
#' @details
#' Matching filenames must contain the patient ID followed by
#' `"Autofluorescence spectra"` and end in `.csv`.
#'
#' @examples
#' \dontrun{
#' latest_af_file(
#'   patient_id = "PT160",
#'   table_spectra_folder = "Output/table_spectra"
#' )
#' }

latest_af_file <- function(patient_id, table_spectra_folder) {
  af_files <- list.files(
    path = table_spectra_folder,
    pattern = paste0(
      patient_id,
      ".*Autofluorescence spectra.*\\.csv$"
    ),
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(af_files) == 0) {
    return(character(0))
  }

  file_info <- file.info(af_files)

  af_files[which.max(file_info$mtime)]
}
