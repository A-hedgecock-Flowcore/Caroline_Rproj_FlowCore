#' Extract a gated flowFrame from a GatingSet
#'
#' Retrieves data from a specified population within one sample of a
#' `GatingSet`, converts the resulting cytoframe to a `flowFrame`, and
#' optionally retains only marker channels.
#'
#' @param gating_set A `flowWorkspace::GatingSet` containing the sample.
#' @param sample_index Numeric index identifying the sample within
#'   `gating_set`.
#' @param population Character string specifying the population to extract.
#'   Defaults to `"singletsSSC"`.
#' @param inverse_transform Logical value indicating whether to apply the
#'   inverse transformation when retrieving the data. Defaults to `FALSE`.
#' @param markers_only Logical value indicating whether to retain only marker
#'   channels. Defaults to `TRUE`.
#'
#' @return A `flowCore::flowFrame` containing the requested gated data.
#'   If `markers_only = TRUE`, the returned object contains only marker
#'   channels.
#'
#' @details
#' The requested population must exist in the sample's gating hierarchy.
#' `sample_index` must identify a valid sample in `gating_set`.
#'
#' @examples
#' \dontrun{
#' ff <- get_gated_flowFrame(
#'   gating_set = PT_gs,
#'   sample_index = 1
#' )
#' }
get_gated_flowFrame <- function(
  gating_set,
  sample_index,
  population = "singletsSSC",
  inverse_transform = FALSE,
  markers_only = TRUE
) {
  gated_cytoframe <- flowWorkspace::gh_pop_get_data(
    gating_set[[sample_index]],
    y = population,
    inverse.transform = inverse_transform
  )

  flow_frame <- flowWorkspace::cytoframe_to_flowFrame(
    gated_cytoframe
  )

  if (markers_only) {
    flow_frame <- flow_frame[, names(markernames(flow_frame))]
  }

  return(flow_frame)
}
