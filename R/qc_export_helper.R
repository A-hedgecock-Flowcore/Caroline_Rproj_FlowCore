#' Make a filesystem-safe filename
#'
#' @param x Character string to sanitize.
#'
#' @return A filesystem-safe character string.
safe_qc_filename <- function(x) {
  x <- as.character(x)

  x <- tools::file_path_sans_ext(
    basename(x)
  )

  x <- gsub(
    "[^A-Za-z0-9_.-]+",
    "_",
    x
  )

  x <- gsub(
    "_+",
    "_",
    x
  )

  x <- gsub(
    "^_|_$",
    "",
    x
  )

  x
}


#' Save a plot-like object as a PNG file
#'
#' Supports ggplot2, patchwork, grob, and other objects that can be printed
#' to a graphics device.
#'
#' @param plot_object Plot-like object to save.
#' @param filename Output PNG filename.
#' @param width Width in inches.
#' @param height Height in inches.
#' @param dpi Resolution in dots per inch.
#'
#' @return Invisibly returns the output filename.
save_qc_plot <- function(
  plot_object,
  filename,
  width = 12,
  height = 8,
  dpi = 160
) {
  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )

  tryCatch(
    {
      ggplot2::ggsave(
        filename = filename,
        plot = plot_object,
        width = width,
        height = height,
        units = "in",
        dpi = dpi,
        bg = "white"
      )
    },
    error = function(ggsave_error) {
      grDevices::png(
        filename = filename,
        width = width,
        height = height,
        units = "in",
        res = dpi
      )

      on.exit(
        grDevices::dev.off(),
        add = TRUE
      )

      tryCatch(
        {
          if (inherits(plot_object, "grob")) {
            grid::grid.newpage()
            grid::grid.draw(plot_object)
          } else {
            print(plot_object)
          }
        },
        error = function(print_error) {
          grDevices::dev.off()

          stop(
            "Could not save QC plot. ggsave error: ",
            conditionMessage(ggsave_error),
            "; print error: ",
            conditionMessage(print_error)
          )
        }
      )
    }
  )

  invisible(filename)
}


#' Save a plot or plot collection as PNG files
#'
#' A single plot-like object is saved as one image. Named lists of plots are
#' saved as one image per list element.
#'
#' @param plot_object Plot-like object or named list of plot-like objects.
#' @param output_folder Directory in which to save images.
#' @param prefix Filename prefix.
#' @param width Width in inches.
#' @param height Height in inches.
#' @param dpi Resolution in dots per inch.
#'
#' @return A character vector of created image paths.
save_qc_plot_collection <- function(
  plot_object,
  output_folder,
  prefix,
  width = 12,
  height = 8,
  dpi = 160
) {
  dir.create(
    output_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (is.null(plot_object)) {
    return(character())
  }

  is_named_plot_collection <- is.list(plot_object) &&
    !inherits(
      plot_object,
      c(
        "gg",
        "ggplot",
        "patchwork",
        "grob",
        "gtable"
      )
    )

  if (!is_named_plot_collection) {
    plot_object <- list(plot_object)
  }

  plot_names <- names(plot_object)

  if (is.null(plot_names) || any(!nzchar(plot_names))) {
    plot_names <- paste0(
      "plot_",
      seq_along(plot_object)
    )
  }

  image_paths <- character(
    length(plot_object)
  )

  for (i in seq_along(plot_object)) {
    filename <- paste0(
      safe_qc_filename(prefix),
      "_",
      safe_qc_filename(plot_names[[i]]),
      ".png"
    )

    image_paths[[i]] <- save_qc_plot(
      plot_object = plot_object[[i]],
      filename = file.path(
        output_folder,
        filename
      ),
      width = width,
      height = height,
      dpi = dpi
    )
  }

  image_paths
}
