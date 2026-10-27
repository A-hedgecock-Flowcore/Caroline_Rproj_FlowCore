# Load openCyto gatingTemplate
cleanup_gates <- data.table::fread(CleanUpGatingTemplateFile)
Template <- gatingTemplate(cleanup_gates)

# Create GatingSets
AutoSpectral_unmixed_files <- list.files(
  path = unmixed_folder,
  pattern = "AutoSpectral\\.fcs$",
  full.names = TRUE,
  ignore.case = TRUE
)

unmixed_controls <- list.files(
  path = unmixed_folder,
  pattern = "FastPoisson\\.fcs$",
  full.names = TRUE,
  ignore.case = TRUE
)

## Make sure there are files where they are supposed to be
if (length(AutoSpectral_unmixed_files) == 0) {
  stop("No AutoSpectral-unmixed files found in: ", unmixed_folder)
} else {
  PT_cs <- load_cytoset_from_fcs(
    AutoSpectral_unmixed_files,
    truncate_max_range = FALSE,
    transformation = FALSE
  )
}

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
}

PT_gs <- GatingSet(PT_cs)
controls_gs <- GatingSet(controls_cs)

# Transform data
## PT_gs
parameters <- colnames(PT_gs)
flurophores_only <- parameters[!str_detect(parameters, "Time|FSC|SSC|Index")]

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

## Controls_gs
if (!is.null(controls_gs)) {
  parameters <- colnames(controls_gs)
  flurophores_only <- parameters[!str_detect(parameters, "Time|FSC|SSC|Index")]

  ### transformation
  Biexponential <- flowjo_biexp_trans(
    channelRange = 4096,
    maxValue = 262144,
    pos = 4.5,
    neg = 0,
    widthBasis = -1000
  )
  BiexTransformList <- transformerList(flurophores_only, Biexponential)

  transform(controls_gs, BiexTransformList)
}

# Ungated NxN
PT_ungated_nxn <- Utility_NbyNPlots(
  PT_gs,
  sample.name = "TUBENAME",
  outpath = output_folder,
  ycolumn = "PE-Cy7-A", # CD45
  bins = 500,
  returntype = "patchwork"
)
if (!is.null(controls_gs)) {
  controls_ungated_nxn <- Utility_NbyNPlots(
    controls_gs,
    sample.name = "TUBENAME",
    outpath = output_folder,
    ycolumn = "PE-Cy7-A", # CD45
    bins = 500,
    returntype = "patchwork"
  )
}

# Cleaup Gates
gt_gating(Template, PT_gs)
PT_cleanup_GatingPlots <- Utility_GatingPlots(
  PT_gs,
  sample.name = "TUBENAME",
  gtFile = cleanup_gates,
  returnType = "patchwork",
  removestrings = ".fcs"
)

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

# Save GatingSet
PT_gs_save_path <- file.path(output_folder, "PT_GatingSet")
save_gs(PT_gs, path = PT_gs_save_path)

if (!is.null(controls_gs)) {
  controls_gs_save_path <- file.path(output_folder, "controls_GatingSet")
  save_gs(controls_gs, path = controls_gs_save_path)
}

# NxN after gating
PT_cleanup_nxn_GatingPlots <- Utility_NbyNPlots(
  PT_gs,
  sample.name = "TUBENAME",
  gatesubset = "singletsSSC",
  ycolumn = "PE-Cy7-A",
  bins = 500
)

if (!is.null(controls_gs)) {
  controls_cleanup_nxn_GatingPlots <- Utility_NbyNPlots(
    controls_gs,
    sample.name = "TUBENAME",
    gatesubset = "singletsSSC",
    ycolumn = "PE-Cy7-A",
    bins = 500
  )
}

qc_gating_plots_file <- file.path(
  output_folder,
  "qc_gating_plots.rds"
)

saveRDS(
  list(
    PT_ungated_NxN = PT_ungated_nxn,
    controls_ungated_NxN = controls_ungated_nxn,
    PT_cleanup_gating = PT_cleanup_GatingPlots,
    controls_cleanup_gating = controls_cleanup_GatingPlots,
    PT_cleaned_NxN = PT_cleanup_nxn_GatingPlots,
    controls_cleaned_NxN = controls_cleanup_nxn_GatingPlots
  ),
  file = qc_gating_plots_file
)
