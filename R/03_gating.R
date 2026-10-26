# Load openCyto gatingTemplate
cleanup_gates <- data.table::fread(CleanUpGatingTemplateFile)
Template <- gatingTemplate(cleanup_gates)

# Create GatingSets
AutoSpectral_unmixed_files <- list.files(
  path = unmixed_folder,
  pattern = "AutoSpectral",
  full.names = TRUE
)

unmixed_controls <- list.files(
  path = unmixed_folder,
  pattern = "FastPoisson",
  full.names = TRUE
)

PT_cs <- load_cytoset_from_fcs(
  AutoSpectral_unmixed_files,
  truncate_max_range = FALSE,
  transformation = FALSE
)

controls_cs <- load_cytoset_from_fcs(
  unmixed_controls,
  truncate_max_range = FALSE,
  transformation = FALSE
)

PT_gs <- GatingSet(PT_cs)
controls_gs <- GatingSet(controls_cs)

# Transform data
## PT_gs
parameters <- colnames(PT_gs)
flurophores_only <- parameters[!str_detect(parameters, "Time|FSC|SSC|Index")]
flurophores_only

# transformation
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
parameters <- colnames(controls_gs)
flurophores_only <- parameters[!str_detect(parameters, "Time|FSC|SSC|Index")]
flurophores_only

# transformation
Biexponential <- flowjo_biexp_trans(
  channelRange = 4096,
  maxValue = 262144,
  pos = 4.5,
  neg = 0,
  widthBasis = -1000
)
BiexTransformList <- transformerList(flurophores_only, Biexponential)

transform(controls_gs, BiexTransformList)

# Ungated NxN
ungated_nxn <- Luciernaga::Utility_NbyNPlots(
  PT_gs,
  sample.name = "TUBENAME",
  outpath = output_folder,
  ycolumn = "PE-Cy7-A", # CD45
  bins = 500,
  returntype = "patchwork"
)
