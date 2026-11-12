# Clustering and Dimenstion Reduction steps after QC of data has been verified

# Setup
library(flowWorkspace)
library(Spectre)
library(tidyverse)
library(data.table)

## file.paths
analysisInfo_folder <- file.path("AnalysisInfo")
parameters_file <- list.files(
  analysisInfo_folder,
  "\\.yml$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(parameters_file) != 1) {
  stop(
    "Expected exactly one (1) .yml (YAML) parameter file in ",
    analysisInfo_folder,
    ". Found ",
    length(parameters_file),
    "."
  )
}

analysis_parameters <- yaml::read_yaml(parameters_file)

output_folder <- file.path("Output")
analysis_outputs <- file.path(output_folder, "analysis_outputs")
if (!dir.exists(analysis_outputs)) {
  dir.create(analysis_outputs)
}

spectre_output_folder <- file.path(output_folder, "Spectre_outputs")
dir.create(spectre_output_folder)

# choose latest gatingSet and load
output_folder_contents <- list.dirs(
  output_folder,
  full.names = TRUE,
  recursive = FALSE
)
PT_gatingSets <- output_folder_contents[str_detect(
  output_folder_contents,
  pattern = "PT_GatingSet"
)]

if (length(PT_gatingSets) > 1) {
  PT_gs_path <- sort(PT_gatingSets)[length(PT_gatingSets)]
  cli::cli_alert_success(paste0("Using GatingSet ", basename(PT_gs_path)))
} else {
  PT_gs_path <- PT_gatingSets
}

PT_gs <- load_gs(PT_gs_path)

# Transfer gs to a data.table
cleaned_data <- gs_pop_get_data(PT_gs, "singlets", inverse.transform = FALSE)
samples_names <- sampleNames(PT_gs)

cleaned_data_list <- map(samples_names, \(sample_id) {
  pop_data <- Biobase::exprs(cleaned_data[[sample_id]])
  dt <- as.data.table(pop_data)
  dt[, Sample := sample_id]
  return(dt)
})

cleaned_data <- rbindlist(cleaned_data_list, use.names = TRUE, fill = TRUE)
cleaned_data[, Sample := str_remove(Sample, " AutoSpectral\\.fcs$")]

message(
  "Cleaned data from samples ",
  paste(unique(cleaned_data$Sample), collapse = ", "),
  "."
)

## clean col names to reflect antigens
markers <- markernames(PT_gs)
fluorophores <- intersect(names(markers), names(cleaned_data))
setnames(
  cleaned_data,
  old = fluorophores,
  new = unname(markers[fluorophores])
)

# select relevant columns
marker_cols <- names(cleaned_data)[
  !str_detect(colnames(cleaned_data), "Time|SSC|FSC|Autofluorescence|ID|Sample")
]
sample_col <- names(cleaned_data)[str_detect(colnames(cleaned_data), "Sample")]
cluster_cols <- names(cleaned_data)[
  !str_detect(
    colnames(cleaned_data),
    "Time|SSC|FSC|Autofluorescence|ID|Sample|CD45|DNA|Viability"
  )
]

# FlowSOM
cleaned_data <- run.flowsom(
  dat = cleaned_data,
  use.cols = cluster_cols,
  xdim = 10,
  ydim = 10,
  meta.k = 15,
  verbose = FALSE
)

# subsample data.table for dimRed
## Cells per sample
cleaned_data |>
  mutate(Sample = as_factor(Sample)) |>
  group_by(Sample) |>
  count()

## downsample to 10K cells per sample
set.seed(12) # reproducibility
subsampled_cleaned_data <- cleaned_data |>
  group_by(Sample) |>
  slice_sample(n = 10000) |>
  ungroup() |>
  as.data.table()

## check subsampling
subsampled_cleaned_data |>
  mutate(Sample = as_factor(Sample)) |>
  group_by(Sample) |>
  count()

# Run UMAP
subsampled_cleaned_data <- run.umap(
  dat = subsampled_cleaned_data,
  use.cols = cluster_cols,
  umap.seed = 42
)

# Visualizations
## heatmap of clusters
flowsom_clusters <- do.aggregate(
  cleaned_data,
  cluster_cols,
  by = "FlowSOM_metacluster"
)
make.pheatmap(
  dat = flowsom_clusters,
  sample.col = "FlowSOM_metacluster",
  plot.cols = cluster_cols,
  file.name = "Heatmap by FlowSOM metaclusters.png",
  path = spectre_output_folder
)

## Expression plots
markers <- cleaned_data |>
  pivot_longer(
    cols = all_of(cluster_cols),
    names_to = "Marker",
    values_to = "Expression"
  )

p <- ggplot(
  data = markers,
  aes(x = Expression, y = FlowSOM_metacluster, fill = FlowSOM_metacluster)
) +
  ggridges::geom_density_ridges(
    scale = 0.9,
    rel_min_height = 0.01,
    alpha = 0.8,
    color = "white",
    linewidth = 0.2
  ) +
  facet_wrap(~Marker, scales = "free_x", ncol = 3) +
  labs(
    x = "Marker expression",
    y = "Metacluster",
    title = "Marker distribution by Metacluster",
    fill = "Cluster"
  ) +
  ggridges::theme_ridges() +
  theme(
    axis.text = element_blank(),
  )

ggsave(
  "Marker expression by FlowSOM metacluster.png",
  p,
  width = 10,
  height = 8,
  bg = "white",
  path = spectre_output_folder
)


p <- ggplot(
  data = markers,
  aes(x = Expression, y = Marker, fill = Marker)
) +
  ggridges::geom_density_ridges(
    scale = 0.9,
    rel_min_height = 0.01,
    alpha = 0.8,
    color = "white",
    linewidth = 0.2
  ) +
  facet_wrap(~FlowSOM_metacluster, scales = "free_x", ncol = 5) +
  labs(
    x = "Marker expression",
    y = "Antigen",
    title = "Marker distribution within Metacluster",
    fill = "Antigen"
  ) +
  ggridges::theme_ridges() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none"
  )

ggsave(
  "Marker distribution within FlowSOM metacluster.png",
  p,
  width = 10,
  height = 8,
  bg = "white",
  path = spectre_output_folder
)
