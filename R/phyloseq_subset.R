
#' @title Nest phyloseq objects according to sample variables
#'
#' @description This function tries to mimic \code{tidyr::nest}, subsetting the \code{phyloseq} object based
#' on all combinations of sample data.
#'
#' @param physeq A \linkS4class{phyloseq} object
#' @param variables A character vector containing the names of the variables to nest by, these variables
#' should be discrete (factor, logical, integer, character)
#' @param keep_variables Logical, should the names of the sample variables being used, be kept
#' in the phyloseq objects?
#' @param keep_empty Logical, should sample data combinations applying to no samples be kept? If
#' \code{TRUE}, the coresponding \code{phyloseq} entity will be coded as \code{NULL}.
#' @param extract_list Logical, if \code{TRUE}, only the the \code{phyloseq} list column is returned
#' @param ... Ignored
#'
#' @return A tibble where the variables are in the in the first columns, while the last column is named \code{phyloseq}
#' and contains the corresponding objects. If \code{extract_list=TRUE}, only this last column is returned.
#' \code{phyloseq} and contains the phyloseq object correspodning to the specific combination of sample variables.
#' @details If only one variable is provided, the element of the \code{phyloseq} column is named after this column
#' @export
#' @import tibble
#' @import dplyr
#' @importFrom rlang new_list
#' @examples
#' library(phyloseq)
#' data("soilrep")
#' subdivide_by_environment(soilrep,variables=c('Treatment','warmed','clipped'),keep_variables=FALSE,keep_empty=FALSE)
subdivide_by_environment <- function(physeq,variables,keep_variables=TRUE,keep_empty=TRUE,extract_list=FALSE,...){
if(!is.character(variables)){
  stop('The variables argument must be a character string of the variables to split by')
}
  sample_data <- phyloseq::sample_data(physeq)
if(variables %in% names(sample_data) %>% all() %>% magrittr::not()){
  stop('Some of the given variables does not exist in the phyloseq object')
}
  unique_values <- lapply(X = variables,
                             FUN = function(variable) sample_data[[variable]] %>% unique()
                             )
  names(unique_values) <- variables
  call_environments <- c(unique_values,list(stringsAsFactors = FALSE))
  environments <- do.call(what = expand.grid, args = call_environments)
  tibble_environments <- tibble::as_tibble(environments)
  res <- add_column(tibble_environments, phyloseq = rlang::new_list(nrow(tibble_environments)))
  sample_names <- phyloseq::sample_names(physeq = physeq)
  sample_variables <- phyloseq::sample_variables(physeq = physeq)
  for(i in seq_len(nrow(tibble_environments))){
  row = tibble_environments[i,]
  names(row) <- names(tibble_environments)
  matching_indecies <- lapply(X=names(row), FUN = function(variable_name) {
    which((sample_data[[variable_name]] %in% row[[variable_name]]))
  })
  selected_samples <- Reduce(f = base::intersect,x = matching_indecies)
  if(length(selected_samples) == 0){
    next()
  }
  selected_samples_names <- sample_names[selected_samples]
  subsetted_phyloseq <- phyloseq::prune_samples(samples = selected_samples_names,x = physeq)
  if(!keep_variables){
    sample_data(subsetted_phyloseq) <- sample_variables %in% variables %>% magrittr::not() %>%
      which() %>% dplyr::select(phyloseq::sample_data(subsetted_phyloseq),.)
  }
  res$phyloseq[[i]] <- subsetted_phyloseq
  }
  if(!keep_empty){
res %<>% dplyr::filter(vapply(res$phyloseq,
                              FUN = function(x) is.null(x) %>% magrittr::not(),FUN.VALUE = logical(length = 1)))
  }
  if(length(variables) == 1){
    names(res$phyloseq) = tibble_environments[[1]]
  }
  if(extract_list){
    return(res$phyloseq)
  }else{
  return(res)
  }
}
