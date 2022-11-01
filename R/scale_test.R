#' Scaling Test (A, B, C, or D)
#'
#' @param data A data frame that contains the generation column and trait column.
#' @param
#' @param generation A vector that contains the generations.
#' @param trait A numeric vector that contains the trait value.
#' @return A data frame that summarizes the number of observations, mean, mean
#' of variance, model.
#' @examples
#' scale_test(data, test="A", generation, trait)
#' scale_test(data, test="B", generation, trait)
#' scale_test(data, test="C", generation, trait)
#' scale_test(data, test="D", generation, trait)
scale_test <- function(data, generation, trait){
  result <- data %>%
    group_by(.data[[generation]]) %>%
    summarize(n = n() , mean = mean(.data[[trait]]), variance = var(.data[[trait]])) %>%
    mutate(mean_var = variance/n, weight_var = 1/mean_var, order = c(5,6,3,4,1,2)) %>% # ensure gen is in order
    arrange(order) %>%
    select(!!generation, n ,mean, variance, mean_var, weight_var)

  res <- as.matrix(result[,-1])
  rownames(res) <- c("P1", "P2", "F1", "F2", "BC1", "BC2")

  ############## A scaling test  ##################
  A <- (res["BC1","mean"] * 2) - res["F1","mean"] - res["P1", "mean"]


  Va <- ( 4 * (res["BC1","mean_var"] ^ 2) +
            (res["F1","mean_var"] ^ 2) +
            (res["P1","mean_var"] ^ 2) )

  # Standard error of A
  Sea <- sqrt(Va)    # Standard Error of A


  ta <- A/Sea     # T-test A

  dfA <- res["BC1","n"] +
    res["P1","n"] +
    res["F1","n"]  - 3   # degree of freedom of A
  ttabA <- qt(0.05, dfA, lower.tail = F )

  sta <-if (abs(ta) > ttabA){"*"}else{"ns"}

  if (abs(ta) > ttabA){
    ta1 <- paste("Dosen't support the additive dominance model")
  }else{
    ta1 <- paste("Supports the additive dominance model")
  }

  ############## B scale test  ##################
  B <- (res["BC2","mean"] * 2) - res["F1","mean"] - res["P2", "mean"]


  Vb <- ( 4 * (res["BC2","mean_var"] ^ 2) +
            (res["F1","mean_var"] ^ 2) +
            (res["P2","mean_var"] ^ 2) )

  # Standard error of B
  Seb <- sqrt(Vb)    # Standard Error of B


  tb <- B/Seb     # T-test B

  dfB <- res["BC2","n"] +
    res["P2","n"] +
    res["F1","n"]  - 3   # degree of freedom of B
  ttabB <- qt(0.05, dfB, lower.tail = F )

  stb <-if (abs(tb) > ttabB){"*"}else{"ns"}

  if (abs(tb) > ttabB){
    tb1 <<- paste("Dosen't support the additive dominance model")
  }else{
    tb1 <- paste("Supports the additive dominance model")
  }

  ###################### C scale test  ##################
  C <- (res["F2","mean"] * 4) - (res["F1","mean"] * 2) -
    res["P1", "mean"] - res["P2", "mean"]


  Vc <- ( 16 * (res["F2","mean_var"] ^ 2)) +
    (4 * (res["F1","mean_var"] ^ 2)) +
    (res["P2", "mean_var"] ^ 2) +
    (res["P1", "mean_var"] ^ 2)

  # Standard error of C
  Sec <- sqrt(Vc)    # Standard Error of C


  tc <- C/Sec     # T-test C

  dfC <- res["F2","n"] +
    res["P1","n"] +
    res["P2", "n"] +
    res["F1","n"]  - 4   # degree of freedom of C
  ttabC <- qt(0.05, dfC, lower.tail = F )

  stc <-if (abs(tc) > ttabC){"*"}else{"ns"}

  if (abs(tc) > ttabC){
    tc1 <<- paste("Dosen't support the additive dominance model")
  }else{
    tc1 <- paste("Supports the additive dominance model")
  }
  ####################### D scale test  ##################
  D <- (res["F2","mean"] * 2) - res["BC1","mean"] - res["BC2", "mean"]


  Vd <- ( 4 * (res["F2","mean_var"] ^ 2) +
            4 * (res["BC1","mean_var"] ^ 2) +
            (res["BC2","mean_var"] ^ 2) )

  # Standard error of D
  Sed <- sqrt(Vd)    # Standard Error of D


  td <- D/Sed     # T-test D

  dfD <- res["F2","n"] +
    res["BC1","n"] +
    res["BC2","n"]  - 3   # degree of freedom of D

  ttabD <- qt(0.05, dfD, lower.tail = F )

  std <-if (abs(td) > ttabD){"*"}else{"ns"}

  if (abs(td) > ttabD){
    td1 <- paste("Dosen't support the additive dominance model")
  }else{
    td1 <- paste("Supports the additive dominance model")
  }
  y <-  tibble(
    "Scaling Test" = c("A", "B", "C", "D"),
    T_Value = c(round(ta,4), round(tb, 4), round(tc, 4), round(td, 4)),
    Significance = c(sta, stb, stc, std),
    comment = c(ta1, tb1, tc1, td1)
  )


  return(y)


}
