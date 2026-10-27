# Setup
library(CompensAID)

## flowframe function and compensAID
source(file.path("R", "get_gated_flowFrame.R"))
source(file.path("R", "run_compensaid.R"))

# Load patient GatingSet if needed
if (!exists("PT_gs", envir = .GlobalEnv)) {
  PT_gs_save_path <- file.path(
    output_folder,
    "PT_GatingSet"
  )

  PT_gs <- load_gs(PT_gs_save_path)
}

# Load controls GatingSet only if needed and available
if (
  !exists("controls_gs", envir = .GlobalEnv) ||
    is.null(controls_gs)
) {
  controls_gs_save_path <- file.path(
    output_folder,
    "controls_GatingSet"
  )

  if (dir.exists(controls_gs_save_path)) {
    controls_gs <- load_gs(controls_gs_save_path)
  } else {
    controls_gs <- NULL
  }
}
# get FlowFrames
## Extract gated patient flowFrames
PT_flowframes <- map(
  seq_along(PT_gs),
  \(i) {
    get_gated_flowFrame(
      gating_set = PT_gs,
      sample_index = i,
      population = "singletsSSC",
      markers_only = TRUE
    )
  }
)

names(PT_flowframes) <- sampleNames(PT_gs)

# Run compensAID
PT_compensaid_results <- purrr::imap(
  PT_flowframes,
  \(flow_frame, sample_name) {
    run_compensaid(
      flow_frame = flow_frame,
      sample_name = sample_name,
      print_plots = FALSE
    )
  }
)

## Run compensAID for controls only when controls are available
if (!is.null(controls_gs)) {
  control_flowframes <- purrr::map(
    seq_along(controls_gs),
    \(i) {
      get_gated_flowFrame(
        gating_set = controls_gs,
        sample_index = i,
        population = "singletsSSC",
        markers_only = TRUE
      )
    }
  )

  names(control_flowframes) <- sampleNames(controls_gs)

  controls_compensaid_results <- purrr::imap(
    control_flowframes,
    \(flow_frame, sample_name) {
      run_compensaid(
        flow_frame = flow_frame,
        sample_name = sample_name,
        print_plots = FALSE
      )
    }
  )
} else {
  control_flowframes <- list()
  controls_compensaid_results <- list()
}
# Save Results
compensaid_output_folder <- file.path(output_folder, "compensAID")
dir.create(compensaid_output_folder, showWarnings = FALSE, recursive = TRUE)
saveRDS(
  PT_compensaid_results,
  file = file.path(compensaid_output_folder, "PT_compensaid_results.rds")
)
saveRDS(
  controls_compensaid_results,
  file = file.path(compensaid_output_folder, "controls_compensaid_results.rds")
)
