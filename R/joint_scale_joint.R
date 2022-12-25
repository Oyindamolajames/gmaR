#' Joint Scaling Test
#'
#' @param data A data frame that contains the generation column and trait column.
#' @param
#' @param generation A vector that contains the generations.
#' @param trait A numeric vector that contains the trait value.
#' @return A list model, genetic parameter and chi square value
#' @examples
#' joint_scale_test(data, generation, trait)

joint_scale_test <- function(data, generation, trait){
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


  j<-matrix(c( sum((1*1*wvarp1), (1*1*wvarp2) , (1*1*wvarf1),
                   (1*1*wvarf2), (1*1*wvarbc1), (1*1*wvarbc2)),
               sum((1*1*wvarp1), (1*-1*wvarp2) , (1*0*wvarf1),
                   (1*0*wvarf2), (1*0.5*wvarbc1), (1*-0.5*wvarbc2)),
               sum((1*0*wvarp1), (1*0*wvarp2) , (1*1*wvarf1),
                   (1*0.5*wvarf2), (1*0.5*wvarbc1), (1*0.5*wvarbc2)),

               sum((1*1*wvarp1), (-1*1*wvarp2) , (0*1*wvarf1),
                   (0*1*wvarf2), (0.5*1*wvarbc1), (-0.5*1*wvarbc2)),
               sum((1*1*wvarp1), (-1*-1*wvarp2) , (0*0*wvarf1),
                   (0*0*wvarf2), (0.5*0.5*wvarbc1), (-0.5*-0.5*wvarbc2)),
               sum((1*0*wvarp1), (-1*0*wvarp2) , (0*1*wvarf1),
                   (0*0.5*wvarf2),(0.5*0.5*wvarbc1), (-0.5*0.5*wvarbc2)),

               sum((0*1*wvarp1), (0*1*wvarp2) , (1*1*wvarf1),
                   (0.5*1*wvarf2), (0.5*1*wvarbc1), (0.5*1*wvarbc2)),
               sum((0*1*wvarp1), (0*-1*wvarp2) , (1*0*wvarf1),
                   (0.5*0*wvarf2), (0.5*0.5*wvarbc1), (0.5*-0.5*wvarbc2)),
               sum((0*0*wvarp1), (0*0*wvarp2) , (1*1*wvarf1),
                   (0.5*0.5*wvarf2),(0.5*0.5*wvarbc1), (0.5*0.5*wvarbc2))
  ),nrow=3,ncol=3,byrow=T)
  colnames(j)<- c("m","d","h")

  s <- matrix(c(sum(1*mp1*wvarp1, 1*mp2*wvarp2, 1*mf1*wvarf1,
                    1*mf2*wvarf2, 1*mbc1*wvarbc1, 1*mbc2*wvarbc2),

                sum(1*mp1*wvarp1, -1*mp2*wvarp2, 0*mf1*wvarf1,
                    0*mf2*wvarf2, 0.5*mbc1*wvarbc1, -0.5*mbc2*wvarbc2),

                sum(0*mp1*wvarp1, 0*mp2*wvarp2, 1*mf1*wvarf1,
                    0.5*mf2*wvarf2, 0.5*mbc1*wvarbc1, 0.5*mbc2*wvarbc2))
              , nrow = 3, ncol = 1, byrow = T)



  ij <- solve(j)

  genp <- ij %*% s


  gm <-  genp[1,1]

  gd <- genp[2,1]

  gh <- genp[3,1]

  sem <- sqrt(ij[1,1])

  sed <- sqrt(ij[2,2])

  seh <- sqrt(ij[3,3])

  ms <- if(gm/sem >1.96){"*"}else{"ns"}
  md <- if(gd/sed >1.96){"*"}else{"ns"}
  mh <- if(gh/seh >1.96){"*"}else{"ns"}


  genpara <- data.frame("Genetic Parameters" = c("m", "d", "h"),
                        "Value" = c(gm,gd,gh),
                        "Standard Error" = c(sem,sed,seh),
                        "sig"=c(ms,md,mh))


  ep1 <- gm + gd
  ep2 <- gm - gd
  ef1 <- gm + gh
  ef2 <- gm + (0.5*gh)
  eb1 <- gm + (0.5*gd) + (0.5*gh)
  eb2 <- gm - (0.5*gd) + (0.5*gh)

  xp1 <- wvarp1*(mp1-ep1)^2
  xp2 <- wvarp2*(mp2-ep2)^2
  xf1 <- wvarf1*(mf1-ef1)^2
  xf2 <- wvarf2*(mf2-ef2)^2
  xbc1 <- wvarbc1*(mbc1-eb1)^2
  xbc2 <- wvarbc2*(mbc2-eb2)^2

  jst <- sum(c(xp1,xp2,xf1,xf2,xbc1,xbc2))

  modelh <- data.frame(Generation = c("p1", "p2", "f1", "f2", "bc1", "bc2"),
                       Mean = c(mp1, mp2, mf1, mf2, mbc1, mbc2),
                       Variance = c(res["P1","variance"], res["P2","variance"], res["F1","variance"], res["F2","variance"], res["BC1","variance"], res["BC2","variance"]),
                       Mean_of_variance = c(res["P1","mean_var"], res["P2","mean_var"], res["F1","mean_var"], res["F2","mean_var"], res["BC1","mean_var"], res["BC2","mean_var"]),
                       Weight = c(wvarp1, wvarp2, wvarf1, wvarf2, wvarbc1, wvarbc2),
                       m = c("1","1","1","1","1","1"),
                       d = c("1","-1","0","0","1/2","-1/2"),
                       h = c("0","0","1","1/2","1/2","1/2"),
                       Observed = c(mp1, mp2, mf1, mf2, mbc1, mbc2),
                       Expected = c(ep1, ep2, ef1, ef2, eb1, eb2),
                       xsqur = c(xp1, xp2, xf1, xf2, xbc1, xbc2))

  result <- list(model=tibble(modelh),
                 chi_square = jst,
                 genpara = genpara)

  return(result)


}
