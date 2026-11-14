# Clustering and Dimesion Reduction steps after QC of data has been verified

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

spectre_output_folder <- file.path(output_folder, "Spectre_outputs")
dir.create(spectre_output_folder)

num_of_clusters <- analysis_parameters$metacluster

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

if (length(PT_gatingSets) == 0) {
  stop(
    "No Patient gatingSets found in ",
    output_folder,
    "."
  )
}

gs_timestamps <- str_match(
  basename(PT_gatingSets),
  "^PT_GatingSet(_Adjusted)?_(\\d{8}_\\d{6})$"
)[, 3]

latest_timestamp <- sort(gs_timestamps, decreasing = TRUE)

PT_gs_path <- PT_gatingSets[1]

cli::cli_alert_success(
  paste0("Using GatingSet ", basename(PT_gs_path))
)
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
  xdim = 14,
  ydim = 14,
  meta.k = num_of_clusters,
  max.meta = 20,
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
    title = "Expression distribution by Antigen",
    fill = "Cluster"
  ) +
  ggridges::theme_ridges() +
  theme(
    axis.text = element_blank()
  )

ggsave(
  "Antigen expression across metacluster.png",
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
    x = "Expression",
    y = "Antigen",
    title = "Marker distribution per Metacluster",
    fill = "Antigen"
  ) +
  ggridges::theme_ridges() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none"
  )

ggsave(
  "Marker distribution per metacluster.png",
  p,
  width = 10,
  height = 8,
  bg = "white",
  path = spectre_output_folder
)

cluster_events <- cleaned_data |>
  group_by(FlowSOM_metacluster, Sample) |>
  count()

p <- ggplot(
  cluster_events,
  aes(x = FlowSOM_metacluster, y = n, fill = FlowSOM_metacluster)
) +
  geom_col() +
  facet_wrap(~Sample) +
  labs(
    title = "Events per metacluster per patient",
    x = "Metacluster",
    y = "Number of Events",
    fill = "Metacluster"
  ) +
  theme_minimal()

ggsave(
  "Events per metacluster per patient.png",
  p,
  width = 15,
  path = spectre_output_folder
)

p <- ggplot(
  cluster_events,
  aes(x = FlowSOM_metacluster, y = n, fill = FlowSOM_metacluster)
) +
  geom_col() +
  labs(
    title = "Events per metacluster",
    x = "Metacluster",
    y = "Number of Events",
    fill = "Metacluster"
  ) +
  theme_minimal()

ggsave(
  "Events per metacluster.png",
  p,
  width = 10,
  path = spectre_output_folder
)

## UMAP display
cluster_centers <- subsampled_cleaned_data |>
  group_by(FlowSOM_metacluster) |>
  dplyr::summarise(
    UMAP_X = median(UMAP_X, na.rm = TRUE),
    UMAP_Y = median(UMAP_Y, na.rm = TRUE),
    .groups = "drop"
  )

p <- ggplot(
  subsampled_cleaned_data,
  aes(UMAP_X, UMAP_Y, color = FlowSOM_metacluster)
) +
  geom_hex(bins = 300) +
  geom_label(
    data = cluster_centers,
    aes(
      UMAP_X,
      UMAP_Y,
      label = FlowSOM_metacluster,
      color = FlowSOM_metacluster
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  labs(
    title = "UMAP with FlowSOM overlays"
  )


ggsave(
  "UMAP with FlowSOM overlay.png",
  p,
  path = spectre_output_folder
)

p <- p + facet_wrap(~Sample)

ggsave("UMAP per patient.png", p, width = 10, path = spectre_output_folder)

markers <- subsampled_cleaned_data |>
  pivot_longer(
    cols = all_of(cluster_cols),
    names_to = "Marker",
    values_to = "Expression"
  )

p <- ggplot(markers, aes(UMAP_X, UMAP_Y)) +
  stat_summary_hex(
    aes(z = Expression, fill = after_stat(value)),
    bins = 300,
    fun = median
  ) +
  facet_wrap(~Marker) +
  scale_fill_viridis_c(name = "Median Expression") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(size = rel(1.5))
  ) +
  labs(title = "Marker Expression")

ggsave(
  "UMAP Marker Expression.png",
  p,
  path = spectre_output_folder
)

# Summary Tables

cell_counts <- cleaned_data |>
  group_by(Sample) |>
  count() |>
  as.data.table()

table <- create.sumtable(
  dat = cleaned_data,
  sample.col = sample_col,
  pop.col = "FlowSOM_metacluster",
  use.cols = cluster_cols,
  counts = cell_counts
)

write_csv(
  table,
  file = file.path(spectre_output_folder, "Summary Table.csv"),
  col_names = TRUE
)

write_csv(
  subsampled_cleaned_data,
  file = file.path(spectre_output_folder, "subsampled_UMAP.csv"),
  col_names = TRUE
)

write_csv(
  cleaned_data,
  file = file.path(spectre_output_folder, "all_cells_noUMAP.csv"),
  col_names = TRUE
)

withr::with_dir(
  spectre_output_folder,
  write.files(
    cleaned_data,
    "Metacluster",
    divide.by = sample_col,
    write.csv = FALSE,
    write.fcs = TRUE
  )
)

withr::with_dir(
  spectre_output_folder,
  write.files(
    cleaned_data,
    "Metacluster",
    write.csv = FALSE,
    write.fcs = TRUE
  )
)
withr::with_dir(
  spectre_output_folder,
  write.files(
    subsampled_cleaned_data,
    "UMAP",
    divide.by = sample_col,
    write.csv = FALSE,
    write.fcs = TRUE
  )
)
withr::with_dir(
  spectre_output_folder,
  write.files(
    subsampled_cleaned_data,
    "UMAP",
    write.csv = FALSE,
    write.fcs = TRUE
  )
)
