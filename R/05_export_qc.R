# Required Output Dirs
export_required_output_dirs <- c(
  file.path(output_folder, "qc_images"),
  file.path(output_folder, "qc_images", "gating"),
  file.path(output_folder, "qc_images", "gating", "patient"),
  file.path(output_folder, "qc_images", "gating", "controls"),
  file.path(output_folder, "qc_images", "compensAID"),
  file.path(output_folder, "qc_images", "compensAID", "patient"),
  file.path(output_folder, "qc_images", "compensAID", "controls"),
  file.path(output_folder, "qc_tables"),
  file.path(output_folder, "variants"),
  file.path(output_folder, "figure_autofluorescence")
)

purrr::walk(
  export_required_output_dirs,
  dir.create,
  showWarnings = FALSE,
  recursive = TRUE
)

# Helper functions
source(file.path("R", "qc_export_helper.R"))

# load objects if needed

read_rds_if_available <- function(path, label) {
  if (!file.exists(path)) {
    warning(
      label,
      " was not found. Skipping: ",
      path
    )

    return(NULL)
  }

  message(
    "Loading ",
    label,
    " from: ",
    path
  )

  readRDS(path)
}

qc_gating_plots_file <- file.path(
  output_folder,
  "qc_gating_plots.rds"
)

patient_compensaid_file <- file.path(
  output_folder,
  "compensAID",
  "PT_compensaid_results.rds"
)

controls_compensaid_file <- file.path(
  output_folder,
  "compensAID",
  "controls_compensaid_results.rds"
)

# Reuse gating plots if they already exist; otherwise load them.
if (
  !exists("qc_gating_plots", envir = .GlobalEnv) ||
    is.null(qc_gating_plots)
) {
  qc_gating_plots <- read_rds_if_available(
    path = qc_gating_plots_file,
    label = "gating plot results"
  )
}

# Reuse patient compensAID results if they already exist; otherwise load them.
if (
  !exists("PT_compensaid_results", envir = .GlobalEnv) ||
    is.null(PT_compensaid_results)
) {
  PT_compensaid_results <- read_rds_if_available(
    path = patient_compensaid_file,
    label = "patient compensAID results"
  )
}

# Reuse control compensAID results if they already exist; otherwise load them.
if (
  !exists("controls_compensaid_results", envir = .GlobalEnv) ||
    is.null(controls_compensaid_results)
) {
  controls_compensaid_results <- read_rds_if_available(
    path = controls_compensaid_file,
    label = "control compensAID results"
  )
}

# Save gating images
tryCatch(
  {
    gating_patient_image_folder <- file.path(
      output_folder,
      "qc_images",
      "gating",
      "patient"
    )

    gating_control_image_folder <- file.path(
      output_folder,
      "qc_images",
      "gating",
      "controls"
    )

    if (is.null(qc_gating_plots)) {
      warning(
        "Gating plot results are unavailable. ",
        "Skipping gating image export."
      )
    } else {
      save_qc_plot_collection(
        plot_object = qc_gating_plots$PT_ungated_NxN,
        output_folder = gating_patient_image_folder,
        prefix = "patient_ungated_NxN"
      )

      save_qc_plot_collection(
        plot_object = qc_gating_plots$PT_cleanup_gating,
        output_folder = gating_patient_image_folder,
        prefix = "patient_cleanup_gating"
      )

      save_qc_plot_collection(
        plot_object = qc_gating_plots$PT_cleaned_NxN,
        output_folder = gating_patient_image_folder,
        prefix = "patient_cleaned_NxN"
      )

      if (!is.null(qc_gating_plots$controls_ungated_NxN)) {
        save_qc_plot_collection(
          plot_object = qc_gating_plots$controls_ungated_NxN,
          output_folder = gating_control_image_folder,
          prefix = "controls_ungated_NxN"
        )
      }

      if (!is.null(qc_gating_plots$controls_cleanup_gating)) {
        save_qc_plot_collection(
          plot_object = qc_gating_plots$controls_cleanup_gating,
          output_folder = gating_control_image_folder,
          prefix = "controls_cleanup_gating"
        )
      }

      if (!is.null(qc_gating_plots$controls_cleaned_NxN)) {
        save_qc_plot_collection(
          plot_object = qc_gating_plots$controls_cleaned_NxN,
          output_folder = gating_control_image_folder,
          prefix = "controls_cleaned_NxN"
        )
      }
    }
  },
  error = function(e) {
    warning(
      "Gating image export failed: ",
      conditionMessage(e)
    )
  }
)

# Save compensAID images
compensaid_patient_image_folder <- file.path(
  output_folder,
  "qc_images",
  "compensAID",
  "patient"
)

compensaid_control_image_folder <- file.path(
  output_folder,
  "qc_images",
  "compensAID",
  "controls"
)

compensaid_table_folder <- file.path(
  output_folder,
  "qc_tables"
)

dir.create(
  compensaid_patient_image_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  compensaid_control_image_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  compensaid_table_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

export_compensaid_qc <- function(
  results,
  image_folder,
  summary_file
) {
  if (is.null(results) || length(results) == 0) {
    empty_summary <- tibble::tibble(
      sample = character(),
      flagged_pairs = integer(),
      review_status = character()
    )

    readr::write_csv(
      empty_summary,
      summary_file
    )

    return(invisible(empty_summary))
  }

  summary_rows <- purrr::imap_dfr(
    results,
    function(result, sample_name) {
      sample_label <- if (!is.null(result$sample_name)) {
        result$sample_name
      } else {
        sample_name
      }

      flagged_pairs <- result$flagged_pairs

      n_flagged <- if (is.null(flagged_pairs)) {
        0L
      } else {
        nrow(flagged_pairs)
      }

      # Save SSI matrix
      if (!is.null(result$compensaid_matrix)) {
        save_qc_plot(
          plot_object = result$compensaid_matrix,
          filename = file.path(
            image_folder,
            paste0(
              safe_qc_filename(sample_label),
              "_SSI_matrix.png"
            )
          ),
          width = 10,
          height = 8,
          dpi = 160
        )
      }

      # Save flagged-pair diagnostic plots
      if (
        !is.null(result$compensaid_plots) &&
          length(result$compensaid_plots) > 0
      ) {
        purrr::iwalk(
          result$compensaid_plots,
          function(plot_object, plot_index) {
            pair_label <- if (
              !is.null(flagged_pairs) &&
                nrow(flagged_pairs) >= plot_index
            ) {
              paste(
                flagged_pairs$primary.marker[[plot_index]],
                flagged_pairs$secondary.marker[[plot_index]],
                sep = "_"
              )
            } else {
              paste0("pair_", plot_index)
            }

            save_qc_plot(
              plot_object = plot_object,
              filename = file.path(
                image_folder,
                paste0(
                  safe_qc_filename(sample_label),
                  "_",
                  safe_qc_filename(pair_label),
                  ".png"
                )
              ),
              width = 10,
              height = 8,
              dpi = 160
            )
          }
        )
      }

      tibble::tibble(
        sample = sample_label,
        flagged_pairs = n_flagged,
        review_status = if (n_flagged > 0) {
          "Review"
        } else {
          "No flagged pairs"
        }
      )
    }
  )

  readr::write_csv(
    summary_rows,
    summary_file
  )

  summary_rows
}
tryCatch(
  {
    patient_compensaid_summary <- export_compensaid_qc(
      results = PT_compensaid_results,
      image_folder = compensaid_patient_image_folder,
      summary_file = file.path(
        compensaid_table_folder,
        "patient_compensaid_summary.csv"
      )
    )

    controls_compensaid_summary <- export_compensaid_qc(
      results = controls_compensaid_results,
      image_folder = compensaid_control_image_folder,
      summary_file = file.path(
        compensaid_table_folder,
        "controls_compensaid_summary.csv"
      )
    )
  },
  error = function(e) {
    compensaid_export_error <<- e

    warning(
      "compensAID image/table export failed: ",
      conditionMessage(e)
    )
  }
)
# Record status

qc_export_status_file <- file.path(
  output_folder,
  "qc_tables",
  "qc_export_status.csv"
)

readr::write_csv(
  tibble::tibble(
    stage = "05_export_qc",
    status = "completed",
    reason = NA_character_
  ),
  qc_export_status_file
)

message(
  "QC export completed."
)

message(
  "QC images saved under: ",
  file.path(output_folder, "qc_images")
)

message(
  "QC tables saved under: ",
  file.path(output_folder, "qc_tables")
)
