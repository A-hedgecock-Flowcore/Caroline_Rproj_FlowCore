# Load openCyto gatingTemplate
cleanup_gates <- data.table::fread(CleanUpGatingTemplateFile)
Template <- gatingTemplate(cleanup_gates)

# PT
## Create GatingSets
AutoSpectral_unmixed_files <- list.files(
  path = unmixed_folder,
  pattern = "AutoSpectral\\.fcs$",
  full.names = TRUE,
  ignore.case = TRUE
)

## Make sure there are files where they are supposed to be
if (length(AutoSpectral_unmixed_files) == 0) {
  stop(
    "No AutoSpectral-unmixed files found in: ",
    unmixed_folder
  )
}

PT_cs <- load_cytoset_from_fcs(
  AutoSpectral_unmixed_files,
  truncate_max_range = FALSE,
  transformation = FALSE
)

PT_gs <- GatingSet(PT_cs)

# Transform data
## PT_gs
parameters <- colnames(PT_gs)
flurophores_only <- parameters[!str_detect(parameters, "Time|FSC|SSC|AF")]

### transformation
Biexponential <- flowjo_biexp_trans(
  channelRange = 4096,
  maxValue = 262144,
  pos = 4.5,
  neg = 0,
  widthBasis = -1000
)
BiexTransformList <- transformerList(flurophores_only, Biexponential)

transform(PT_gs, BiexTransformList)

# Clean
source("R/get_gated_flowFrame.R")

PT_flowframes <- map(
  seq_along(PT_gs),
  \(i) {
    get_gated_flowFrame(
      gating_set = PT_gs,
      sample_index = i,
      population = "root",
      markers_only = FALSE
    )
  }
)

names(PT_flowframes) <- sampleNames(PT_gs)

PT_flowframes <- purrr::imap(
  PT_flowframes,
  \(flow_frame, sample_name) {
    PeacoQC::RemoveMargins(
      ff = flow_frame,
      channels = parameters,
      output = "frame"
    )
  }
)


PT_PeacoQC_results <- purrr::imap(
  PT_flowframes,
  \(flow_frame, sample_name) {
    PeacoQC::PeacoQC(
      ff = flow_frame,
      channels = flurophores_only,
      determine_good_cells = "all",
      plot = TRUE,
      save_fcs = FALSE,
      output_directory = output_folder
    )
  }
)

clean_ffs <- map(PT_PeacoQC_results, "FinalFF")
names(clean_ffs) <- names(PT_PeacoQC_results)

PT_fs <- flowSet(frames = clean_ffs)
PT_gs <- GatingSet(PT_fs)

# Ungated NxN
PT_ungated_nxn <- Utility_NbyNPlots(
  PT_gs,
  sample.name = "TUBENAME",
  outpath = output_folder,
  ycolumn = "SSC-A",
  returntype = "patchwork"
)

# Cleaup Gates
gt_gating(Template, PT_gs)

PT_cleanup_GatingPlots <- Utility_GatingPlots(
  PT_gs,
  sample.name = "TUBENAME",
  gtFile = cleanup_gates,
  returnType = "patchwork",
  removestrings = ".fcs"
)

# NxN after gating
PT_cleanup_nxn_GatingPlots <- Utility_NbyNPlots(
  PT_gs,
  sample.name = "TUBENAME",
  gatesubset = "singlets",
  ycolumn = "SSC-A"
)


# Save GatingSet
GatingSetTimestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
PT_gs_save_path <- file.path(
  output_folder,
  paste0("PT_GatingSet", "_", GatingSetTimestamp)
)
save_gs(PT_gs, path = PT_gs_save_path)


# Controls
## GatingSet
unmixed_controls <- list.files(
  path = unmixed_folder,
  pattern = "FastPoisson\\.fcs$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(unmixed_controls) == 0) {
  warning(
    "No FastPoisson control files found in: ",
    unmixed_folder,
    ". Control gating and control compensAID will be skipped."
  )

  controls_gs <- NULL
} else {
  controls_cs <- load_cytoset_from_fcs(
    unmixed_controls,
    truncate_max_range = FALSE,
    transformation = FALSE
  )

  controls_gs <- GatingSet(controls_cs)
}

## Controls_gs
if (!is.null(controls_gs)) {
  parameters <- colnames(controls_gs)
  flurophores_only <- parameters[
    !str_detect(parameters, "Time|FSC|SSC|Index")
  ]

  ### transformation
  Biexponential <- flowjo_biexp_trans(
    channelRange = 4096,
    maxValue = 262144,
    pos = 4.5,
    neg = 0,
    widthBasis = -500
  )
  BiexTransformList <- transformerList(flurophores_only, Biexponential)

  transform(controls_gs, BiexTransformList)
}

# Ungated NxN
if (!is.null(controls_gs)) {
  controls_ungated_nxn <- Utility_NbyNPlots(
    controls_gs,
    sample.name = "TUBENAME",
    outpath = output_folder,
    ycolumn = "SSC-A",
    returntype = "patchwork"
  )
}

# Gates
if (!is.null(controls_gs)) {
  gt_gating(Template, controls_gs)

  controls_cleanup_GatingPlots <- Utility_GatingPlots(
    controls_gs,
    sample.name = "TUBENAME",
    gtFile = cleanup_gates,
    returnType = "patchwork",
    removestrings = ".fcs"
  )
}

# NxN after gating

if (!is.null(controls_gs)) {
  controls_cleanup_nxn_GatingPlots <- Utility_NbyNPlots(
    controls_gs,
    sample.name = "TUBENAME",
    gatesubset = "singlets",
    ycolumn = "SSC-A"
  )
}

# Save gatingset

if (!is.null(controls_gs)) {
  controls_gs_save_path <- file.path(
    output_folder,
    paste0(
      "controls_GatingSet",
      "_",
      GatingSetTimestamp
    )
  )
  save_gs(controls_gs, path = controls_gs_save_path)
}

# save plots
qc_gating_plots_file <- file.path(
  output_folder,
  "qc_gating_plots.rds"
)

qc_gating_plots <- list(
  PT_ungated_NxN = PT_ungated_nxn,
  PT_cleanup_gating = PT_cleanup_GatingPlots,
  PT_cleaned_NxN = PT_cleanup_nxn_GatingPlots
)

if (!is.null(controls_gs)) {
  qc_gating_plots$controls_ungated_NxN <- controls_ungated_nxn
  qc_gating_plots$controls_cleanup_gating <- controls_cleanup_GatingPlots
  qc_gating_plots$controls_cleaned_NxN <- controls_cleanup_nxn_GatingPlots
}

qc_gating_plots_file <- file.path(
  output_folder,
  "qc_gating_plots.rds"
)

saveRDS(
  qc_gating_plots,
  file = qc_gating_plots_file
)
