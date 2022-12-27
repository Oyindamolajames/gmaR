#' Computes Genetic Parameters
#'
#' @param data A data frame that contains the generation column and trait column.
#' @param generation A vector that contains the generations.
#' @param trait A numeric vector that contains the trait value.
#' @return A data frame that summarizes heritabilty, heterosis, pcv, gcv, etc.
#' @examples
#' get_summary(data, generation, trait)


genetic_parameters <- function(data, generation, trait){
  result <- data %>%
    group_by(.data[[generation]]) %>%
    summarize(n = n() , mean = mean(.data[[trait]]), variance = var(.data[[trait]])) %>%
    mutate(mean_var = variance/n, weight_var = 1/mean_var, order = c(5,6,3,4,1,2)) %>% # ensure gen is in order
    arrange(order) %>%
    select(!!generation, n ,mean, variance, mean_var, weight_var)

  res <- as.matrix(result[,-1])
  rownames(res) <- c("P1", "P2", "F1", "F2", "BC1", "BC2")

  wvarp1 <- res["P1","weight_var"]
  wvarp2 <- res["P2","weight_var"]
  wvarf1 <- res["F1","weight_var"]
  wvarf2 <- res["F2","weight_var"]
  wvarbc1 <- res["BC1","weight_var"]
  wvarbc2 <- res["BC2","weight_var"]

  mp1 <- res["P1","mean"]
  mp2 <- res["P2","mean"]
  mf1 <- res["F1","mean"]
  mf2 <- res["F2","mean"]
  mbc1 <- res["BC1","mean"]
  mbc2 <- res["BC2","mean"]

  vp1 <- res["P1","variance"]
  vp2 <- res["P2","variance"]
  vf1 <- res["F1","variance"]
  vf2 <- res["F2","variance"]
  vbc1 <- res["BC1","variance"]
  vbc2 <- res["BC2","variance"]

  meanvarp1 <- res["P1","mean_var"]
  meanvarp2 <- res["P2","mean_var"]
  meanvarf1 <- res["F1","mean_var"]
  meanvarf2 <- res["F2","mean_var"]
  meanvarbc1 <- res["BC1","mean_var"]
  meanvarbc2 <- res["BC2","mean_var"]

  mean_trait <- data %>%
    summarize(mean = mean(.data[[trait]])) %>%
    pull(mean)


  #Heritability
  VE <- (vp1 + vp2 + (vf1 * 2)) / 4
  VG <- vf2 - VE
  VP <- VE + VG
  VD <- 4 * (vbc1 + vbc2 - vf2 - VE)
  VA <- 2 * (vf2 - (0.25 * VD) - VE)
  HB <- ((0.5 * VA) + (0.25 * VD) / (0.5 * VA) + (0.25 * VD) + VE)
  HN <- ((0.5 * VA) / (0.5 * VA) + (0.25 * VD) + VE)
  ht <- mf1-((mp1+mp2)/2)
  pcv <- sqrt(VG+VE)/mean_trait *100
  gcv <- sqrt(VG)/mean_trait *100
  ga<-round((VG/VP)*2.06*sqrt(VP),4)
  gam<-round((ga/mean_trait) *100,4)
  mph <- mf1-((mp1+mp2)/2)
  bph <- if (mp1>mp2){
    mf1-mp1
  }else{
    mf1-mp2
  }
  id <- mf1-mf2

  heri<- data.frame(Components= c("VP","VE","VG","VA","VD","HB","HN","GA","GAM","heterosis",
                                  "PCV","GCV","MPH","BPH","ID"),
                    Value=c(vf2,VE,VG,VA,VD,round(HB,2),round(HN,2),ga,gam,ht,pcv,gcv,mph,bph,id))

  return(heri)
}
