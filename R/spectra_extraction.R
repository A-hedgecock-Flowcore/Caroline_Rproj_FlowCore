#' Spectra Extraction
#'
#' Searches the AutoSpectral output directory for an existing automated
#' spectra CSV file. If multiple matching files exist, the most recently
#' modified file is loaded. If no matching file exists, new spectra are
#' generated using `AutoSpectral::get.spectra.automated()`.
#'
#' @param output_folder Character string specifying the main output directory.
#' @param control.dir Character string specifying the AutoSpectral control
#'   directory.
#' @param control_filename Character string specifying the AutoSpectral
#'   control definition file.
#' @param asp AutoSpectral object required by `get.spectra.automated()`.
#'
#' @return A spectra matrix returned by `AutoSpectral::read.spectra()` or
#'   `get.spectra.automated()`.
#'
#' @details
#' Existing files must match the pattern:
#' `Automated_autospectral_spectra_YYYYMMDD_HHMMSS.csv`.
#'
#' @examples
#' \dontrun{
#' spectra <- load_or_create_spectra(
#'   output_folder = "Output",
#'   control.dir = "Input/SCC",
#'   control_filename = "Input/AnalysisInfo/fcs_control_file",
#'   asp = asp
#' )
#' }
spectra_extraction <- function(
  output_folder,
  control_dir,
  control_filename,
  asp
) {
  if (base::length(existing_spectra_files) > 0) {
    file_info <- file.info(existing_spectra_files)

    latest_file <- rownames(
      file_info[which.max(file_info$mtime), , drop = FALSE]
    )

    spectra <- read.spectra(
      spectra.file = basename(latest_file),
      spectra.dir = dirname(latest_file)
    )
    cli::cli_alert_success(paste0(
      "Using spectra file: ",
      basename(latest_file),
      " ."
    ))
  } else {
    # if no file exists, create
    full_control_dir <- normalizePath(control.dir, mustWork = TRUE)
    full_control_filename <- normalizePath(control_filename, mustWork = TRUE)

    spectra <- withr::with_dir(
      output_folder,
      get.spectra.automated(
        control.dir = full_control_dir,
        control.def.file = full_control_filename,
        asp = asp,
        figures = TRUE,
        verbose = FALSE
      )
    )
    cli::cli_alert_success("Created new spectra.")
  }
  return(spectra)
}
