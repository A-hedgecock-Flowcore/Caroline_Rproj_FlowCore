# AnalysisPipeline.R
#
# Top-level driver for:
#   01_setup.R
#   02_unmixing.R
#   03_gating.R
#   04_compensAID.R
#
# Run from the project root:
#   source("AnalysisPipeline.R")
#
# Or from a terminal:
#   Rscript AnalysisPipeline.R

options(
  stringsAsFactors = FALSE,
  warn = 1
)

pipeline_start_time <- Sys.time()

pipeline_start_timestamp <- format(
  pipeline_start_time,
  "%Y%m%d_%H%M%S"
)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!dir.exists(file.path(project_root, "AnalysisInfo"))) {
  stop(
    "This script must be run from the project root. ",
    "Could not find AnalysisInfo/ in: ",
    project_root
  )
}

if (file.exists(file.path(project_root, "renv", "activate.R"))) {
  source(file.path(project_root, "renv", "activate.R"))
}

required_directories <- c(
  "AnalysisInfo",
  "Input",
  "Input/SCC",
  "Input/Samples",
  "R",
  "Output"
)

missing_directories <- required_directories[
  !dir.exists(file.path(project_root, required_directories))
]

if (length(missing_directories) > 0) {
  stop(
    "The following required directories are missing: ",
    paste(missing_directories, collapse = ", ")
  )
}

log_directory <- file.path(project_root, "Output", "logs")

dir.create(
  log_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

pipeline_log_file <- file.path(
  log_directory,
  paste0(
    "pipeline_",
    pipeline_start_timestamp,
    ".log"
  )
)

pipeline_log_connection <- file(
  pipeline_log_file,
  open = "wt"
)

sink(pipeline_log_connection, type = "output")
sink(pipeline_log_connection, type = "message")

pipeline_error <- NULL

tryCatch(
  {
    message("Starting analysis pipeline")
    message("Project root: ", project_root)
    message("Start time: ", pipeline_start_time)

    message("Running setup")
    source(file.path(project_root, "R", "01_setup.R"))

    message(
      "Running Autospectral unmixing with per-patient, per-cell AF extraction and fluorophore variants"
    )
    source(file.path(project_root, "R", "02_unmixing.R"))

    message("Running automated cleanup gating")
    source(file.path(project_root, "R", "03_gating.R"))

    message("Running compensAID package for unmixing evaluation")
    source(file.path(project_root, "R", "04_compensAID.R"))

    message("Pipeline stages completed")
  },
  error = function(e) {
    pipeline_error <<- e

    message(
      "PIPELINE ERROR: ",
      conditionMessage(e)
    )

    message(
      "Review the log file: ",
      pipeline_log_file
    )
  },
  finally = {
    sink(type = "message")
    sink(type = "output")
    close(pipeline_log_connection)
  }
)

pipeline_end_time <- Sys.time()

pipeline_status <- if (is.null(pipeline_error)) {
  "completed"
} else {
  "failed"
}

pipeline_summary <- list(
  status = pipeline_status,
  start_time = pipeline_start_time,
  end_time = pipeline_end_time,
  duration_minutes = as.numeric(
    difftime(
      pipeline_end_time,
      pipeline_start_time,
      units = "mins"
    )
  ),
  project_root = project_root,
  log_file = pipeline_log_file,
  error = if (is.null(pipeline_error)) {
    NULL
  } else {
    list(
      class = class(pipeline_error),
      message = conditionMessage(pipeline_error)
    )
  },
  input = list(
    Patient_samples = if (exists("Patient_samples")) {
      Patient_samples
    } else {
      character()
    },
    Patient_unstained_samples = if (exists("Patient_unstained_samples")) {
      Patient_unstained_samples
    } else {
      character()
    },
    nonpatient_samples = if (exists("nonpatient_samples")) {
      nonpatient_samples
    } else {
      character()
    }
  ),
  output = list(
    output_folder = if (exists("output_folder")) {
      output_folder
    } else {
      file.path(project_root, "Output")
    },
    unmixed_folder = if (exists("unmixed_folder")) {
      unmixed_folder
    } else {
      file.path(project_root, "Output", "Unmixed_files")
    }
  ),
  status_files = list(
    patient_unmixing = if (exists("patient_unmix_status_file")) {
      patient_unmix_status_file
    } else {
      NA_character_
    },
    nonpatient_unmixing = if (exists("nonpatient_unmix_status_file")) {
      nonpatient_unmix_status_file
    } else {
      NA_character_
    },
    af_status = if (exists("af_status_file")) {
      af_status_file
    } else {
      NA_character_
    },
    qc_gating_plots = file.path(
      project_root,
      "Output",
      "qc_gating_plots.rds"
    ),
    patient_compensaid = file.path(
      project_root,
      "Output",
      "compensAID",
      "PT_compensaid_results.rds"
    ),
    controls_compensaid = file.path(
      project_root,
      "Output",
      "compensAID",
      "controls_compensaid_results.rds"
    )
  )
)

summary_output_folder <- if (exists("output_folder")) {
  output_folder
} else {
  file.path(project_root, "Output")
}

dir.create(
  summary_output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  pipeline_summary,
  file = file.path(
    summary_output_folder,
    "pipeline_summary.rds"
  )
)

message(
  "Pipeline status: ",
  pipeline_status
)

message(
  "Pipeline duration: ",
  round(pipeline_summary$duration_minutes, 2),
  " minutes"
)

message(
  "Pipeline summary saved to: ",
  file.path(summary_output_folder, "pipeline_summary.rds")
)

if (!is.null(pipeline_error)) {
  stop(pipeline_error)
}
