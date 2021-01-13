#' @importFrom stats prcomp
#'
#' @title
#' Plot trajectory for time series
#'
#' @description
#' When looking at PCA-plots for different time series, we often want to study their evolution
#' in time. This convenience functions does so by adding arrows between the samples in their
#' chronological order in the PCA plot
#'
#' @param time_series_list A list of \link{OTU_time_series} objects
#' @param distance The name of any distance metric in \link{vegdist}
#' @param subset Character, the subsets of time series to plot, based on the
#' names of the list.
#' @param label Logical, should each point be labelled with the time point?
#' @param label_size Positive numeric, the font size of the label. Ignored if \code{text=FALSE}.
#' @param color Character, the color mapping to use in the plot. Defaults to the time series
#' @param linetype Character, the linetype mapping to be used. Defaults to solid
#' @details The function utilzes the entire table to make the principal components,
#' even if some data points are excluded. Note the columns \code{'time_series'} and \code{'time_points'}
#' which are available in the object being returned.
#' @return A \link{ggplot} object showing the trajectory in the two first principal components
#' @importFrom vegan vegdist
#' @import ggfortify
#' @examples
#' library(micInt)
#' library(phyloseq)
#' data("seawater")
#' phyloseq_list <- subdivide_by_environment(seawater,"Reactor")
#' time_series <- OTU_time_series(phyloseq_list,"Week")
#' plot_trajectory(time_series,distance = "euclidean")
#' @export
plot_trajectory <- function(time_series_list, distance = "bray",
                            subset = names(time_series_list), label = FALSE, label_size = 3,
                            color = NULL, linetype = NULL) {
  if (!is.list(time_series_list) || !all(vapply(
    X = time_series_list,
    function(x) inherits(x, "OTU_time_series"),
    logical(length = 1)
  )
  )
  ) {
    stop("Input must be a list of elements of class OTU_time_series")
  }
  if (time_series_list %>% names() %>% is.null()) {
    stop("The list must be named")
  }
  frames <- lapply(names(time_series_list), function(name) {
    table <- time_series_list[[name]]@table
    time_points <- time_series_list[[name]]@time_points
    # Sorts the observations in chronological order
    table <- table[order(time_points, decreasing = FALSE), ]
    table$time_series <- name
    table$time_points <- sort(time_points)
    return(table)
  })
  names(frames) <- names(time_series_list)
  # Combines the different frames in order to make a PCA ordination
  all_frames <- do.call(rbind, frames)
  dist_matrix <- vegdist(all_frames[, colnames(all_frames) %in% c("time_series", "time_points") %>% `!`()] %>% as.matrix(),
    method = distance, diag = TRUE, upper = TRUE
  )
  PCA_all <- prcomp(dist_matrix)
  fort_table_all <- fortify(PCA_all)
  fort_table_all$time_series <- all_frames$time_series
  fort_table_all$time_points <- all_frames$time_points
  fort_table <- fort_table_all[fort_table_all$time_series %in% subset, ]
  # Finding the proprotion explained by the principal components
  total_variance <- sum(PCA_all$sdev^2)
  PC1_proportion <- PCA_all$sdev[1]^2 / total_variance
  PC2_proportion <- PCA_all$sdev[2]^2 / total_variance
  if (is.null(color)) {
    color <- "time_series"
  }
  if (is.null(linetype)) {
    res <- ggplot(data = fort_table, mapping = aes_string(
      x = "PC1",
      y = "PC2",
      group = "time_series",
      color = color
    )) + geom_point(size = 0)
  }
  else {
    res <- ggplot(data = fort_table, mapping = aes_string(
      x = "PC1",
      y = "PC2",
      group = "time_series",
      color = color,
      linetype = linetype
    )) + geom_point(size = 0)
  }
  if (label) {
    res <- res + geom_text(aes_string(label = "time_points"), size = label_size, color = "black")
  }
  if (is.null(color)) {
    res <- res + guides(color = guide_legend("Time series"))
  }
  res <- res + geom_path(aes_string(linetype = linetype),
    arrow = arrow(length = unit(0.1, "inches"))
  )
  res <- res + xlab(sprintf("PC1 (%.2f%%)", PC1_proportion * 100)) +
    ylab(sprintf("PC2 (%.2f%%)", PC2_proportion * 100))
  return(res)
}
