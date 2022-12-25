#' Perfect Fit Test
#'
#' @param data A data frame that contains the generation column and trait column.
#' @param
#' @param generation A vector that contains the generations.
#' @param trait A numeric vector that contains the trait value.
#' @return A list model, genetic parameter and chi square value
#' @examples
#' perfect_fit(data, generation, trait)

perfectf_fit <- function(data, generation, trait) {
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

  meanvarp1 <- res["P1","mean_var"]
  meanvarp2 <- res["P2","mean_var"]
  meanvarf1 <- res["F1","mean_var"]
  meanvarf2 <- res["F2","mean_var"]
  meanvarbc1 <- res["BC1","mean_var"]
  meanvarbc2 <- res["BC2","mean_var"]




  #Perfect fit

  pm <- 0.5*mp1 + 0.5*mp2 + 4*mf1 -2*mbc1 - 2*mbc2
  pd <- 0.5*mp1 - 0.5*mp2
  ph <- 6*mbc1 + 6*mbc2- 8*mf2 - mf1 -1.5*mp1 - 1.5*mp2
  pi <- 2*mbc1 + 2*mbc2 - 4*mf2
  pj <- 2*mbc1 - mf1 - 2*mbc2 + mp2
  pl <- mp1 + mp2 + 2*mf1 + 4*mf2 - 4*mbc1 - 4*mbc2

  #se_pm <- sqrt(0.25*meanvarp1 + 0.25*mp2 + 16*meanvarf1 +4*meanvarbc1 + 4*meanvarbc2)
  #se_pd <-sqrt(0.25*meanvarp1 + 0.25*meanvarp2)
  #se_ph <- sqrt(36*meanvarbc1 + 36*meanvarbc2+ 64*meanvarf2 + meanvarf1 + 1.5*meanvarp1 + 2.25*meanvarp2)
  #se_pi <- sqrt(4*meanvarbc1 + 4*meanvarbc2 + 16*meanvarf2)
  #se_pj <- sqrt(4*meanvarbc1 + meanvarf1 + 4*meanvarbc2 + meanvarp2)
  #se_pl <-sqrt(meanvarp1 + meanvarp2 + 4*meanvarf1 + 16*meanvarf2 + 16*meanvarbc1 + 16*meanvarbc2)

  se_pm <- sqrt(0.25*(meanvarp1^2) + 0.25*(mp2^2) + 16*(meanvarf1^2) +4*(meanvarbc1^2) + 4*(meanvarbc2^2))
  se_pd <-sqrt(0.25*(meanvarp1^2) + 0.25*(meanvarp2^2))
  se_ph <- sqrt(36*(meanvarbc1^2) + 36*(meanvarbc2^2)+ 64*(meanvarf2^2) + (meanvarf1^2) + 1.5*(meanvarp1^2) + 2.25*(meanvarp2^2))
  se_pi <- sqrt(4*(meanvarbc1^2) + 4*(meanvarbc2^2) + 16*(meanvarf2^2))
  se_pj <- sqrt(4*(meanvarbc1^2) + (meanvarf1^2) + 4*(meanvarbc2^2) + (meanvarp2^2))
  se_pl <-sqrt((meanvarp1^2) + (meanvarp2^2) + 4*(meanvarf1^2) + 16*(meanvarf2^2) + 16*(meanvarbc1^2) + 16*(meanvarbc2^2))




  tpm <- pm/se_pm
  tpd <- pd/se_pd
  tph <- ph/se_ph
  tpi <- pi/se_pi
  tpj <- pj/se_pj
  tpl <- pl/se_pl

  p_m <- if(abs(tpm)>1.96){
    pm
  }else{0}

  p_d <- if(abs(tpd)>1.96){
    pd
  }else{0}

  p_h <- if(abs(tph)>1.96){
    ph
  }else{0}

  p_i <- if(abs(tpi)>1.96){
    pi
  }else{0}

  p_j <- if(abs(tpj)>1.96){
    pj
  }else{0}

  p_l <- if(abs(tpl)>1.96){
    pl
  }else{0}
  #significant
  p_ms <- if(abs(tpm)>1.96){
    "*"
  }else{"ns"}

  p_ds <- if(abs(tpd)>1.96){
    "*"
  }else{"ns"}

  p_hs <- if(abs(tph)>1.96){
    "*"
  }else{"ns"}

  p_is <- if(abs(tpi)>1.96){
    "*"
  }else{"ns"}

  p_js <- if(abs(tpj)>1.96){
    "*"
  }else{"ns"}

  p_ls <- if(abs(tpl)>1.96){
    "*"
  }else{"ns"}


  pp1 <- p_m + p_d + p_i
  pp2 <- p_m - p_d + p_i
  pf1 <- p_m + p_h + p_l
  pf2 <- p_m + 0.5*p_h + 0.25*p_l
  pbc1 <- p_m + 0.5*p_d +0.5* p_h + 0.25*p_i + 0.25*p_j + 0.25*p_l
  pbc2<-  p_m - 0.5*p_d +0.5* p_h + 0.25*p_i - 0.25*p_j + 0.25*p_l


  pxp1 <- wvarp1*(mp1-pp1)^2
  pxp2 <- wvarp2*(mp2-pp2)^2
  pxf1 <- wvarf1*(mf1-pf1)^2
  pxf2 <- wvarf2*(mf2-pf2)^2
  pxbc1 <- wvarbc1*(mbc1-pbc1)^2
  pxbc2 <- wvarbc2*(mbc2-pbc2)^2

  px <- sum(c(pxp1,pxp2,pxf1,pxf2,pxbc1,pxbc2))

  pf <- "Perfect fit "

  pfg <- data.frame(GEN =c("m","d","h","i","j","l"),
                    VALUE =c(pm,pd,ph,pi,pj,pl),
                    Stand_err= c(se_pm,se_pd,se_ph,se_pi,se_pj,se_pl),
                    Val= c(p_m,p_d,p_h,p_i,p_j,p_l),
                    Sig= c(p_ms,p_ds,p_hs,p_is,p_js,p_ls)
  )
  pfh<- "perfect fit model"
  pfm <- data.frame(Gen= c("p1","p2","f1","f2","bc1","bc2"),
                    Mean = c(mp1, mp2, mf1, mf2, mbc1, mbc2),
                    expected=c(pp1,pp2,pf1,pf2,pbc1,pbc2),
                    chisqr=c(pxp1,pxp2,pxf1,pxf2,pxbc1,pxbc2))
  pfx <- c("The Chiquare value is",px)

  result <- list("Genetic Parameters" = pfg,
                 "model" = pfm,
                 "chi_square = pfx")

  return (result)
}
