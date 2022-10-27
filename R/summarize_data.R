#' Summarize data
#'
#' @param data A data frame that contains the generation column and trait column.
#' @param generation A vector that contains the generations.
#' @param trait A numeric vector that contains the trait value.
#' @return A data frame that summarizes the number of observations, mean, mean
#' of variance, model.
#' @examples
#' get_summary(data, generation, trait)

get_summary <- function(data, generation, trait){
  result <- data %>%
    group_by(.data[[generation]]) %>%
    summarize(n = n() , mean = mean(.data[[trait]]), variance = var(.data[[trait]])) %>%
    mutate(mean_var = variance/n, weight_var = 1/mean_var, order = c(5,6,3,4,1,2)) %>%
    arrange(order) %>%
    select(n, !!generation, mean, variance, mean_var, weight_var)

  return(result)
}
