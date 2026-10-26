#' Run compensAID quality control on a flowFrame
#'
#' Runs compensAID on one gated flowFrame, generates the SSI matrix,
#' identifies channel pairs with SSI values below a specified threshold, and
#' creates diagnostic plots for flagged channel pairs.
#'
#' @param flow_frame A `flowCore::flowFrame` containing the gated sample data.
#' @param sample_name Optional character string identifying the sample.
#'   Defaults to `"Unnamed sample"`.
#' @param ssi_threshold Numeric threshold used to identify flagged channel
#'   pairs. Defaults to `-1`.
#' @param print_plots Logical value indicating whether flagged-pair plots
#'   should be printed immediately. Defaults to `TRUE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{sample_name}{The sample name used during processing.}
#'   \item{compensaid_results}{The results returned by
#'   `CompensAID::CompensAID()`.}
#'   \item{compensaid_matrix}{The SSI matrix plot returned by
#'   `CompensAID::PlotMatrix()`.}
#'   \item{flagged_pairs}{A tibble listing channel and marker pairs with SSI
#'   values below `ssi_threshold`.}
#'   \item{compensaid_plots}{A list of diagnostic plots for flagged pairs.}
#' }
#'
#' @details
#' Marker names are matched to flow cytometry parameter names using the
#' `markernames` and `colnames` metadata stored in `flow_frame`.
#'
#' If no SSI values fall below `ssi_threshold`, `flagged_pairs` is returned as
#' an empty tibble and `compensaid_plots` is returned as an empty list.
#'
#' @examples
#' \dontrun{
#' result <- run_compensaid(
#'   flow_frame = PT_flowframes[[1]],
#'   sample_name = names(PT_flowframes)[1]
#' )
#'
#' result$flagged_pairs
#' }
run_compensaid <- function(
  flow_frame,
  sample_name = NULL,
  ssi_threshold = -1,
  print_plots = FALSE
) {
  if (is.null(sample_name)) {
    sample_name <- "Unnamed sample"
  }

  message("Running compensAID for: ", sample_name)

  compensaid_results <- CompensAID::CompensAID(flow_frame)

  compensaid_matrix <- CompensAID::PlotMatrix(
    compensaid_results
  )

  ssi_matrix <- compensaid_results$matrix

  flagged_index <- which(
    ssi_matrix < ssi_threshold,
    arr.ind = TRUE
  )

  if (nrow(flagged_index) == 0) {
    message(
      "No SSI values below ",
      ssi_threshold,
      " found for: ",
      sample_name
    )

    flagged_pairs <- tibble::tibble(
      primary.channel = character(),
      secondary.channel = character(),
      primary.marker = character(),
      secondary.marker = character()
    )

    compensaid_plots <- list()
  } else {
    channel_names <- rownames(ssi_matrix)
    marker_names <- flowCore::markernames(flow_frame)
    parameter_names <- flowCore::colnames(flow_frame)

    flagged_pairs <- tibble::tibble(
      primary.channel = channel_names[
        flagged_index[, "col"]
      ],
      secondary.channel = channel_names[
        flagged_index[, "row"]
      ]
    ) |>
      dplyr::mutate(
        primary.marker = marker_names[
          match(primary.channel, parameter_names)
        ],
        secondary.marker = marker_names[
          match(secondary.channel, parameter_names)
        ]
      )

    compensaid_plots <- purrr::pmap(
      flagged_pairs[
        c("primary.marker", "secondary.marker")
      ],
      \(primary.marker, secondary.marker) {
        CompensAID::PlotDotSSI(
          compensaid_results,
          flow_frame,
          primary = primary.marker,
          secondary = secondary.marker,
          showScores = TRUE
        )
      }
    )

    if (print_plots) {
      purrr::walk(
        compensaid_plots,
        print
      )
    }
  }

  list(
    sample_name = sample_name,
    compensaid_results = compensaid_results,
    compensaid_matrix = compensaid_matrix,
    flagged_pairs = flagged_pairs,
    compensaid_plots = compensaid_plots
  )
}
