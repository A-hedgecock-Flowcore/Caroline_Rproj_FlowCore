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
