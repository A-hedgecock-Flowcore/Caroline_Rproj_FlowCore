# Spectral Cytometry settings
asp <- get.autospectral.param(cytometer = "aurora")
control.dir <- single_color_control_folder
control_filename <- list.files(
  path = output_folder,
  pattern = analysis_parameters$AutoSpectral_control_file_name,
  full.names = TRUE,
  ignore.case = TRUE
)
