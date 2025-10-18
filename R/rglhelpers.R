#' Convert tidy segment data to interleaved point format
#'
#' Takes a data frame / matrix of line segments where each row defines a start and end point
#' (e.g. `x`, `y`, `z`, `xend`, `yend`, `zend`) and returns a long-form, interleaved
#' data frame where each segment contributes two rows — one for the start and one for
#' the end point — preserving any non-coordinate metadata.
#'
#' @param data A data frame containing segment start and end coordinates.
#' @param coord Character vector of coordinate names (default: `c("x", "y", "z")`).
#' @param end_suffix String suffix appended to coordinate names for the end point
#'   (default: `"end"`).
#'
#' @return A data frame with one row per point, containing:
#'   - `point`: `"start"` or `"end"`
#'   - the coordinate columns (`coord`)
#'   - any additional metadata replicated for both points.
#'
#' @examples
#' df <- data.frame(x = 0, y = 0, z = 0, xend = 1, yend = 1, zend = 1, id = 1)
#' convert_tidy_segments_to_interleaved(df)
#'
#' @seealso [geom_segment()] for ggplot2-style segment data.
#' @export
convert_tidy_segments_to_interleaved <- function(data, coord = c("x", "y", "z"), end_suffix = "end") {
  # Columns assumed to be coordinates
  start_cols <- coord
  end_cols <- paste0(coord, end_suffix)

  # Safety check
  missing_cols <- setdiff(c(start_cols, end_cols), colnames(data))
  if (length(missing_cols)) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  n <- nrow(data)

  # Point labels
  point <- rep(c("start", "end"), times = n)

  # Build interleaved coordinates
  xyz_mat <- t(cbind(data[, start_cols, drop=FALSE], data[,end_cols,drop=FALSE]))
  xyz_mat <- matrix(xyz_mat, ncol = length(coord), byrow = TRUE)
  xyz_data <- as.data.frame(xyz_mat, stringsAsFactors = FALSE)
  colnames(xyz_data) <- coord

  #browser()
  # Keep and replicate non-coordinate (meta) columns
  meta_cols <- setdiff(colnames(data), c(start_cols, end_cols))
  meta_data <- if (length(meta_cols)) {
    data[rep(seq_len(n), each = 2), meta_cols, drop = FALSE]
  } else {
    NULL
  }

  # Assemble output (segment info + meta + coords)
  newcoords <- data.frame(
    point = point,
    xyz_data,
    row.names = NULL,
    check.names = FALSE
  )

  if(!is.null(meta_data)){
    newcoords <- cbind(newcoords, meta_data)
  }

  newcoords
}
