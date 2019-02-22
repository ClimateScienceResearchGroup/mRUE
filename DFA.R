library(readr)
library(MARSS)

mRUE_data <- read_csv("D:/desktop/mRUE_data_v2.csv")
#View(mRUE_data)


all_dat <- mRUE_data



Z_est_T1 = array(, dim=c(3,2,8))
lm_r2_T1 = array(, dim=c(1,2,8))
lm_coef_T1 = array(, dim=c(1,2,8))
nan_count_T1 = array(, dim=c(3,2,8))

Z_est_T2 = array(, dim=c(3,2,8))
lm_r2_T2 = array(, dim=c(1,2,8))
lm_coef_T2 = array(, dim=c(1,2,8))
nan_count_T2 = array(, dim=c(3,2,8))

Z_est_T3 = array(, dim=c(3,2,8))
lm_r2_T3 = array(, dim=c(1,2,8))
lm_coef_T3 = array(, dim=c(1,2,8))
nan_count_T3 = array(, dim=c(3,2,8))

Z_est_T4 = array(, dim=c(3,2,8))
lm_r2_T4 = array(, dim=c(1,2,8))
lm_coef_T4 = array(, dim=c(1,2,8))
nan_count_T4 = array(, dim=c(3,2,8))

Z_est_T5 = array(, dim=c(3,2,8))
lm_r2_T5 = array(, dim=c(1,2,8))
lm_coef_T5 = array(, dim=c(1,2,8))
nan_count_T5 = array(, dim=c(3,2,8))

Z_est_T6 = array(, dim=c(3,2,8))
lm_r2_T6 = array(, dim=c(1,2,8))
lm_coef_T6 = array(, dim=c(1,2,8))
nan_count_T6 = array(, dim=c(3,2,8))

Z_est_T7 = array(, dim=c(3,2,8))
lm_r2_T7 = array(, dim=c(1,2,8))
lm_coef_T7 = array(, dim=c(1,2,8))
nan_count_T7 = array(, dim=c(3,2,8))
#########Loop through all years


######## START OF T1 LOOP
for (z in 1:8){
  year<-2008+z
  ## use only the 10 years from 1980-1989
  yr_frst <- year
  yr_last <- year
  mRUE_dat <- all_dat[all_dat[, "Year"] >= yr_frst & all_dat[, 
                                                              "Year"] <= yr_last, ]
  ## create vector of group names
  obs_data <- c("T1_day_WUEs", "T1_day_LUEs", "T1_day_CUEs")
  
  
  ## get only the RUEs
  dat_working <- mRUE_dat[, obs_data]
  NEP <- mRUE_dat[,"T1_day_NEP"]
  
  ## transpose data so time goes across columns
  dat_working <- t(dat_working)
  
  #### NAN counting
  nan_count_T1[1,1,z] = sum(is.na(dat_working[1,]))
  nan_count_T1[2,1,z] = sum(is.na(dat_working[2,]))
  nan_count_T1[3,1,z] = sum(is.na(dat_working[3,]))
  
  ## get number of time series
  N_ts <- dim(dat_working)[1]
  ## get length of time series
  TT <- dim(dat_working)[2]
  
  y_bar <- apply(dat_working, 1, mean, na.rm = TRUE)
  dat <- dat_working - y_bar
  rownames(dat) <- rownames(dat_working)

  clr <- c("blue", "darkgreen", "darkred")

  
  ## number of processes
  mm <- 1
  ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
  BB <- "identity"  # diag(mm)
  ## 'uu' is a column vector of 0's
  uu <- "zero"  # matrix(0,mm,1)
  ## 'CC' and 'cc' are for covariates
  CC <- "zero"  # matrix(0,mm,1)
  cc <- "zero"  # matrix(0,1,wk_last)
  ## 'QQ' is identity
  QQ <- "identity"  # diag(mm)

  
  ## 'ZZ' is loadings matrix, 1x1
  Z_vals <- list("z11", "z21", "z31")
  ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 1, byrow = TRUE)
  ZZ
  
  ## 'aa' is the offset/scaling
  aa <- "zero"
  ## 'DD' and 'd' are for covariates
  DD <- "zero"  # matrix(0,mm,1)
  dd <- "zero"  # matrix(0,1,wk_last)
  ## 'RR' is var-cov matrix for obs errors
  RR <- "diagonal and unequal"
  
  ## list with specifications for model vectors/matrices
  mod_list <- list(B = BB, U = uu, C = CC, c = cc, Q = QQ, Z = ZZ, 
                   A = aa, D = DD, d = dd, R = RR)
  ## list with model inits
  init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
  ## list with model control parameters
  con_list <- list(maxit = 10000, allow.degen = TRUE)
  
  ## fit MARSS
  dfa_1 <- MARSS(y = dat, model = mod_list, inits = init_list, 
                 control = con_list)
  
  ## get the estimated ZZ
  Z_est <- coef(dfa_1, type = "matrix")$Z
  
  ## get the inverse of the rotation matrix
  #H_inv <- varimax(Z_est)$rotmat
  
  ## rotate factor loadings
  #Z_rot = Z_est %*% H_inv
  ## rotate processes
  proc_rot = Z_est %*% dfa_1$states
  
  
  
  ###### two stage loading plots
  ylbl <- obs_data
  w_ts <- seq(dim(dat)[2])
  layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
  ## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
  par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))
  ## plot the processes
  for (i in 1:mm) {
    ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
    ## set up plot area
    plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
         xlab = "", ylab = "", xaxt = "n")
    ## draw zero-line
    abline(h = 0, col = "gray")
    ## plot trend line
    lines(w_ts, proc_rot[i, ], lwd = 2)
    lines(w_ts, proc_rot[i, ], lwd = 2)
    ## add panel labels
    mtext(paste("State", i), side = 3, line = 0.5)
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
  }
  ## plot the loadings
  minZ <- 0
  ylm <- c(-1, 1) * max(abs(Z_est))
  for (i in 1:mm) {
    plot(c(1:N_ts)[abs(Z_est[, i]) > minZ], as.vector(Z_est[abs(Z_est[, 
                                                                      i]) > minZ, i]), type = "h", lwd = 2, xlab = "", ylab = "", 
         xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col = clr)
    for (j in 1:N_ts) {
      if (Z_est[j, i] > minZ) {
        text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
             col = clr[j])
      }
      if (Z_est[j, i] < -minZ) {
        text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
             col = clr[j])
      }
      abline(h = 0, lwd = 1.5, col = "gray")
    }
    mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
  }
  
  

  
  get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
    ## empty list for results
    fits <- list()
    ## extra stuff for var() calcs
    Ey <- MARSS:::MARSShatyt(MLEobj)
    ## model params
    ZZ <- coef(MLEobj, type = "matrix")$Z
    ## number of obs ts
    nn <- dim(Ey$ytT)[1]
    ## number of time steps
    TT <- dim(Ey$ytT)[2]
    ## get the inverse of the rotation matrix
    #H_inv <- varimax(ZZ)$rotmat
    ## check for covars
    if (!is.null(dd)) {
      DD <- coef(MLEobj, type = "matrix")$D
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states + DD %*% dd
    } else {
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states
    }
    ## Var in model fits
    VtT <- MARSSkfss(MLEobj)$VtT
    VV <- NULL
    for (tt in 1:TT) {
      RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
                                                           , tt] %*% t(ZZ)
      SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
        t(MLEobj$states[, tt, drop = FALSE])
      VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
    }
    SE <- sqrt(VV)
    ## upper & lower (1-alpha)% CI
    fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
    fits$lo <- qnorm(alpha/2) * SE + fits$ex
    return(fits)
  }
  
  
  ## get model fits & CI's
  mod_fit <- get_DFA_fits(dfa_1)
  ## plot the fits
  ylbl <- obs_data
  par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
                                                               0, 0, 0))
  for (i in 1:N_ts) {
    up <- mod_fit$up[i, ]
    mn <- mod_fit$ex[i, ]
    lo <- mod_fit$lo[i, ]
    plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
         cex.lab = 1.2, ylim = c(min(lo), max(up)))
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
    lines(w_ts, dat[i, ], pch = 16, col = clr[i], lwd = 2)
    lines(w_ts, up, col = "darkgray")
    lines(w_ts, mn, col = "black", lwd = 2)
    lines(w_ts, lo, col = "darkgray")
    
  }
  
  ############## Regression work
  
  tempNEP <- array(unlist(NEP), dim = c(1, dim(NEP)))
  
  # Fit regression model
  sat.mod1 <- lm(tempNEP[1:365] ~ proc_rot[1,1:365]) # regression formula
  #sat.mod2 <- lm(tempNEP[1:365] ~ proc_rot[2,1:365]) # regression formula
  
  
  ################# SAVING DATA EACH YEAR
  Z_est_T1[,,z] <- Z_est
  lm_r2_T1[1,1,z] <- summary(sat.mod1)$r.squared
  lm_coef_T1[1,1,z] <- summary(sat.mod1)$coefficients[2, 1]
  lm_coef_T1[1,2,z] <- summary(sat.mod1)$coefficients[1, 1]
}
######## END OF T1 LOOP







######## START OF T2 LOOP
for (z in 1:8){
  year<-2008+z
  ## use only the 10 years from 1980-1989
  yr_frst <- year
  yr_last <- year
  mRUE_dat <- all_dat[all_dat[, "Year"] >= yr_frst & all_dat[, 
                                                              "Year"] <= yr_last, ]
  ## create vector of group names
  obs_data <- c("T2_day_WUEs", "T2_day_LUEs", "T2_day_CUEs")
  
  
  ## get only the RUEs
  dat_working <- mRUE_dat[, obs_data]
  NEP <- mRUE_dat[,"T2_day_NEP"]
  
  ## transpose data so time goes across columns
  dat_working <- t(dat_working)
  
  #### NAN counting
  nan_count_T2[1,1,z] = sum(is.na(dat_working[1,]))
  nan_count_T2[2,1,z] = sum(is.na(dat_working[2,]))
  nan_count_T2[3,1,z] = sum(is.na(dat_working[3,]))
  
  ## get number of time series
  N_ts <- dim(dat_working)[1]
  ## get length of time series
  TT <- dim(dat_working)[2]
  
  y_bar <- apply(dat_working, 1, mean, na.rm = TRUE)
  dat <- dat_working - y_bar
  rownames(dat) <- rownames(dat_working)
  

  clr <- c("blue", "darkgreen", "darkred")

  
  ## number of processes
  mm <- 1
  ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
  BB <- "identity"  # diag(mm)
  ## 'uu' is a column vector of 0's
  uu <- "zero"  # matrix(0,mm,1)
  ## 'CC' and 'cc' are for covariates
  CC <- "zero"  # matrix(0,mm,1)
  cc <- "zero"  # matrix(0,1,wk_last)
  ## 'QQ' is identity
  QQ <- "identity"  # diag(mm)
  
  
  ## 'ZZ' is loadings matrix, 1x1
  Z_vals <- list("z11", "z21", "z31")
  ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 1, byrow = TRUE)
  ZZ
  
  ## 'aa' is the offset/scaling
  aa <- "zero"
  ## 'DD' and 'd' are for covariates
  DD <- "zero"  # matrix(0,mm,1)
  dd <- "zero"  # matrix(0,1,wk_last)
  ## 'RR' is var-cov matrix for obs errors
  RR <- "diagonal and unequal"
  
  ## list with specifications for model vectors/matrices
  mod_list <- list(B = BB, U = uu, C = CC, c = cc, Q = QQ, Z = ZZ, 
                   A = aa, D = DD, d = dd, R = RR)
  ## list with model inits
  init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
  ## list with model control parameters
  con_list <- list(maxit = 10000, allow.degen = TRUE)
  
  ## fit MARSS
  dfa_1 <- MARSS(y = dat, model = mod_list, inits = init_list, 
                 control = con_list)
  
  ## get the estimated ZZ
  Z_est <- coef(dfa_1, type = "matrix")$Z
  
  ## get the inverse of the rotation matrix
  #H_inv <- varimax(Z_est)$rotmat
  
  ## rotate factor loadings
  #Z_rot = Z_est %*% H_inv
  ## rotate processes
  proc_rot = Z_est %*% dfa_1$states
  
  
  
  ###### two stage loading plots
  ylbl <- obs_data
  w_ts <- seq(dim(dat)[2])
  layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
  ## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
  par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))
  ## plot the processes
  for (i in 1:mm) {
    ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
    ## set up plot area
    plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
         xlab = "", ylab = "", xaxt = "n")
    ## draw zero-line
    abline(h = 0, col = "gray")
    ## plot trend line
    lines(w_ts, proc_rot[i, ], lwd = 2)
    lines(w_ts, proc_rot[i, ], lwd = 2)
    ## add panel labels
    mtext(paste("State", i), side = 3, line = 0.5)
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
  }
  ## plot the loadings
  minZ <- 0
  ylm <- c(-1, 1) * max(abs(Z_est))
  for (i in 1:mm) {
    plot(c(1:N_ts)[abs(Z_est[, i]) > minZ], as.vector(Z_est[abs(Z_est[, 
                                                                      i]) > minZ, i]), type = "h", lwd = 2, xlab = "", ylab = "", 
         xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col = clr)
    for (j in 1:N_ts) {
      if (Z_est[j, i] > minZ) {
        text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
             col = clr[j])
      }
      if (Z_est[j, i] < -minZ) {
        text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
             col = clr[j])
      }
      abline(h = 0, lwd = 1.5, col = "gray")
    }
    mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
  }
  
  
  
  get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
    ## empty list for results
    fits <- list()
    ## extra stuff for var() calcs
    Ey <- MARSS:::MARSShatyt(MLEobj)
    ## model params
    ZZ <- coef(MLEobj, type = "matrix")$Z
    ## number of obs ts
    nn <- dim(Ey$ytT)[1]
    ## number of time steps
    TT <- dim(Ey$ytT)[2]
    ## get the inverse of the rotation matrix
    #H_inv <- varimax(ZZ)$rotmat
    ## check for covars
    if (!is.null(dd)) {
      DD <- coef(MLEobj, type = "matrix")$D
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states + DD %*% dd
    } else {
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states
    }
    ## Var in model fits
    VtT <- MARSSkfss(MLEobj)$VtT
    VV <- NULL
    for (tt in 1:TT) {
      RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
                                                           , tt] %*% t(ZZ)
      SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
        t(MLEobj$states[, tt, drop = FALSE])
      VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
    }
    SE <- sqrt(VV)
    ## upper & lower (1-alpha)% CI
    fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
    fits$lo <- qnorm(alpha/2) * SE + fits$ex
    return(fits)
  }
  
  
  ## get model fits & CI's
  mod_fit <- get_DFA_fits(dfa_1)
  ## plot the fits
  ylbl <- obs_data
  par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
                                                               0, 0, 0))
  for (i in 1:N_ts) {
    up <- mod_fit$up[i, ]
    mn <- mod_fit$ex[i, ]
    lo <- mod_fit$lo[i, ]
    plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
         cex.lab = 1.2, ylim = c(min(lo), max(up)))
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
    lines(w_ts, dat[i, ], pch = 16, col = clr[i], lwd = 2)
    lines(w_ts, up, col = "darkgray")
    lines(w_ts, mn, col = "black", lwd = 2)
    lines(w_ts, lo, col = "darkgray")
    
  }
  
  ############## Regression work
  
  tempNEP <- array(unlist(NEP), dim = c(1, dim(NEP)))
  
  # Fit regression model
  sat.mod1 <- lm(tempNEP[1:365] ~ proc_rot[1,1:365]) # regression formula
  #sat.mod2 <- lm(tempNEP[1:365] ~ proc_rot[2,1:365]) # regression formula
  
  
  ################# SAVING DATA EACH YEAR
  Z_est_T2[,,z] <- Z_est
  lm_r2_T2[1,1,z] <- summary(sat.mod1)$r.squared
  lm_coef_T2[1,1,z] <- summary(sat.mod1)$coefficients[2, 1]
  lm_coef_T2[1,2,z] <- summary(sat.mod1)$coefficients[1, 1]
}
######## END OF T2 LOOP




######## START OF T3 LOOP
for (z in 1:8){
  year<-2008+z
  ## use only the 10 years from 1980-1989
  yr_frst <- year
  yr_last <- year
  mRUE_dat <- all_dat[all_dat[, "Year"] >= yr_frst & all_dat[, 
                                                              "Year"] <= yr_last, ]
  ## create vector of group names
  obs_data <- c("T3_day_WUEs", "T3_day_LUEs", "T3_day_CUEs")
  
  
  ## get only the RUEs
  dat_working <- mRUE_dat[, obs_data]
  NEP <- mRUE_dat[,"T3_day_NEP"]
  
  ## transpose data so time goes across columns
  dat_working <- t(dat_working)
  
  #### NAN counting
  nan_count_T3[1,1,z] = sum(is.na(dat_working[1,]))
  nan_count_T3[2,1,z] = sum(is.na(dat_working[2,]))
  nan_count_T3[3,1,z] = sum(is.na(dat_working[3,]))
  
  ## get number of time series
  N_ts <- dim(dat_working)[1]
  ## get length of time series
  TT <- dim(dat_working)[2]
  
  y_bar <- apply(dat_working, 1, mean, na.rm = TRUE)
  dat <- dat_working - y_bar
  rownames(dat) <- rownames(dat_working)
  

  clr <- c("blue", "darkgreen", "darkred")

  
  ## number of processes
  mm <- 1
  ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
  BB <- "identity"  # diag(mm)
  ## 'uu' is a column vector of 0's
  uu <- "zero"  # matrix(0,mm,1)
  ## 'CC' and 'cc' are for covariates
  CC <- "zero"  # matrix(0,mm,1)
  cc <- "zero"  # matrix(0,1,wk_last)
  ## 'QQ' is identity
  QQ <- "identity"  # diag(mm)

  Z_vals <- list("z11", "z21", "z31")
  ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 1, byrow = TRUE)
  ZZ
  
  ## 'aa' is the offset/scaling
  aa <- "zero"
  ## 'DD' and 'd' are for covariates
  DD <- "zero"  # matrix(0,mm,1)
  dd <- "zero"  # matrix(0,1,wk_last)
  ## 'RR' is var-cov matrix for obs errors
  RR <- "diagonal and unequal"
  
  ## list with specifications for model vectors/matrices
  mod_list <- list(B = BB, U = uu, C = CC, c = cc, Q = QQ, Z = ZZ, 
                   A = aa, D = DD, d = dd, R = RR)
  ## list with model inits
  init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
  ## list with model control parameters
  con_list <- list(maxit = 10000, allow.degen = TRUE)
  
  ## fit MARSS
  dfa_1 <- MARSS(y = dat, model = mod_list, inits = init_list, 
                 control = con_list)
  
  ## get the estimated ZZ
  Z_est <- coef(dfa_1, type = "matrix")$Z
  
  ## get the inverse of the rotation matrix
  #H_inv <- varimax(Z_est)$rotmat
  
  ## rotate factor loadings
  #Z_rot = Z_est %*% H_inv
  ## rotate processes
  proc_rot = Z_est %*% dfa_1$states
  
  
  
  ###### two stage loading plots
  ylbl <- obs_data
  w_ts <- seq(dim(dat)[2])
  layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
  ## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
  par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))
  ## plot the processes
  for (i in 1:mm) {
    ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
    ## set up plot area
    plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
         xlab = "", ylab = "", xaxt = "n")
    ## draw zero-line
    abline(h = 0, col = "gray")
    ## plot trend line
    lines(w_ts, proc_rot[i, ], lwd = 2)
    lines(w_ts, proc_rot[i, ], lwd = 2)
    ## add panel labels
    mtext(paste("State", i), side = 3, line = 0.5)
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
  }
  ## plot the loadings
  minZ <- 0
  ylm <- c(-1, 1) * max(abs(Z_est))
  for (i in 1:mm) {
    plot(c(1:N_ts)[abs(Z_est[, i]) > minZ], as.vector(Z_est[abs(Z_est[, 
                                                                      i]) > minZ, i]), type = "h", lwd = 2, xlab = "", ylab = "", 
         xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col = clr)
    for (j in 1:N_ts) {
      if (Z_est[j, i] > minZ) {
        text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
             col = clr[j])
      }
      if (Z_est[j, i] < -minZ) {
        text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
             col = clr[j])
      }
      abline(h = 0, lwd = 1.5, col = "gray")
    }
    mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
  }
  
  
  
  get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
    ## empty list for results
    fits <- list()
    ## extra stuff for var() calcs
    Ey <- MARSS:::MARSShatyt(MLEobj)
    ## model params
    ZZ <- coef(MLEobj, type = "matrix")$Z
    ## number of obs ts
    nn <- dim(Ey$ytT)[1]
    ## number of time steps
    TT <- dim(Ey$ytT)[2]
    ## get the inverse of the rotation matrix
    #H_inv <- varimax(ZZ)$rotmat
    ## check for covars
    if (!is.null(dd)) {
      DD <- coef(MLEobj, type = "matrix")$D
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states + DD %*% dd
    } else {
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states
    }
    ## Var in model fits
    VtT <- MARSSkfss(MLEobj)$VtT
    VV <- NULL
    for (tt in 1:TT) {
      RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
                                                           , tt] %*% t(ZZ)
      SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
        t(MLEobj$states[, tt, drop = FALSE])
      VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
    }
    SE <- sqrt(VV)
    ## upper & lower (1-alpha)% CI
    fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
    fits$lo <- qnorm(alpha/2) * SE + fits$ex
    return(fits)
  }
  
  
  ## get model fits & CI's
  mod_fit <- get_DFA_fits(dfa_1)
  ## plot the fits
  ylbl <- obs_data
  par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
                                                               0, 0, 0))
  for (i in 1:N_ts) {
    up <- mod_fit$up[i, ]
    mn <- mod_fit$ex[i, ]
    lo <- mod_fit$lo[i, ]
    plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
         cex.lab = 1.2, ylim = c(min(lo), max(up)))
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
    lines(w_ts, dat[i, ], pch = 16, col = clr[i], lwd = 2)
    lines(w_ts, up, col = "darkgray")
    lines(w_ts, mn, col = "black", lwd = 2)
    lines(w_ts, lo, col = "darkgray")
    
  }
  
  ############## Regression work
  
  tempNEP <- array(unlist(NEP), dim = c(1, dim(NEP)))
  
  # Fit regression model
  sat.mod1 <- lm(tempNEP[1:365] ~ proc_rot[1,1:365]) # regression formula
  #sat.mod2 <- lm(tempNEP[1:365] ~ proc_rot[2,1:365]) # regression formula

  
  ################# SAVING DATA EACH YEAR
  Z_est_T3[,,z] <- Z_est
  lm_r2_T3[1,1,z] <- summary(sat.mod1)$r.squared
  lm_coef_T3[1,1,z] <- summary(sat.mod1)$coefficients[2, 1]
  lm_coef_T3[1,2,z] <- summary(sat.mod1)$coefficients[1, 1]
}
######## END OF T3 LOOP





######## START OF T4 LOOP
for (z in 1:8){
  year<-2008+z
  ## use only the 10 years from 1980-1989
  yr_frst <- year
  yr_last <- year
  mRUE_dat <- all_dat[all_dat[, "Year"] >= yr_frst & all_dat[, 
                                                              "Year"] <= yr_last, ]
  ## create vector of group names
  obs_data <- c("T4_day_WUEs", "T4_day_LUEs", "T4_day_CUEs")
  
  
  ## get only the RUEs
  dat_working <- mRUE_dat[, obs_data]
  NEP <- mRUE_dat[,"T4_day_NEP"]
  
  ## transpose data so time goes across columns
  dat_working <- t(dat_working)
  
  #### NAN counting
  nan_count_T4[1,1,z] = sum(is.na(dat_working[1,]))
  nan_count_T4[2,1,z] = sum(is.na(dat_working[2,]))
  nan_count_T4[3,1,z] = sum(is.na(dat_working[3,]))
  
  ## get number of time series
  N_ts <- dim(dat_working)[1]
  ## get length of time series
  TT <- dim(dat_working)[2]
  
  y_bar <- apply(dat_working, 1, mean, na.rm = TRUE)
  dat <- dat_working - y_bar
  rownames(dat) <- rownames(dat_working)
  
  clr <- c("blue", "darkgreen", "darkred")

  
  ## number of processes
  mm <- 1
  ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
  BB <- "identity"  # diag(mm)
  ## 'uu' is a column vector of 0's
  uu <- "zero"  # matrix(0,mm,1)
  ## 'CC' and 'cc' are for covariates
  CC <- "zero"  # matrix(0,mm,1)
  cc <- "zero"  # matrix(0,1,wk_last)
  ## 'QQ' is identity
  QQ <- "identity"  # diag(mm)
  
  
  ## 'ZZ' is loadings matrix, 1x1
  Z_vals <- list("z11", "z21", "z31")
  ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 1, byrow = TRUE)
  ZZ
  
  ## 'aa' is the offset/scaling
  aa <- "zero"
  ## 'DD' and 'd' are for covariates
  DD <- "zero"  # matrix(0,mm,1)
  dd <- "zero"  # matrix(0,1,wk_last)
  ## 'RR' is var-cov matrix for obs errors
  RR <- "diagonal and unequal"
  
  ## list with specifications for model vectors/matrices
  mod_list <- list(B = BB, U = uu, C = CC, c = cc, Q = QQ, Z = ZZ, 
                   A = aa, D = DD, d = dd, R = RR)
  ## list with model inits
  init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
  ## list with model control parameters
  con_list <- list(maxit = 10000, allow.degen = TRUE)
  
  ## fit MARSS
  dfa_1 <- MARSS(y = dat, model = mod_list, inits = init_list, 
                 control = con_list)
  
  ## get the estimated ZZ
  Z_est <- coef(dfa_1, type = "matrix")$Z
  
  ## get the inverse of the rotation matrix
  #H_inv <- varimax(Z_est)$rotmat
  
  ## rotate factor loadings
  #Z_rot = Z_est %*% H_inv
  ## rotate processes
  proc_rot = Z_est %*% dfa_1$states
  
  
  
  ###### two stage loading plots
  ylbl <- obs_data
  w_ts <- seq(dim(dat)[2])
  layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
  ## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
  par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))
  ## plot the processes
  for (i in 1:mm) {
    ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
    ## set up plot area
    plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
         xlab = "", ylab = "", xaxt = "n")
    ## draw zero-line
    abline(h = 0, col = "gray")
    ## plot trend line
    lines(w_ts, proc_rot[i, ], lwd = 2)
    lines(w_ts, proc_rot[i, ], lwd = 2)
    ## add panel labels
    mtext(paste("State", i), side = 3, line = 0.5)
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
  }
  ## plot the loadings
  minZ <- 0
  ylm <- c(-1, 1) * max(abs(Z_est))
  for (i in 1:mm) {
    plot(c(1:N_ts)[abs(Z_est[, i]) > minZ], as.vector(Z_est[abs(Z_est[, 
                                                                      i]) > minZ, i]), type = "h", lwd = 2, xlab = "", ylab = "", 
         xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col = clr)
    for (j in 1:N_ts) {
      if (Z_est[j, i] > minZ) {
        text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
             col = clr[j])
      }
      if (Z_est[j, i] < -minZ) {
        text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
             col = clr[j])
      }
      abline(h = 0, lwd = 1.5, col = "gray")
    }
    mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
  }
  

  
  get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
    ## empty list for results
    fits <- list()
    ## extra stuff for var() calcs
    Ey <- MARSS:::MARSShatyt(MLEobj)
    ## model params
    ZZ <- coef(MLEobj, type = "matrix")$Z
    ## number of obs ts
    nn <- dim(Ey$ytT)[1]
    ## number of time steps
    TT <- dim(Ey$ytT)[2]
    ## get the inverse of the rotation matrix
    #H_inv <- varimax(ZZ)$rotmat
    ## check for covars
    if (!is.null(dd)) {
      DD <- coef(MLEobj, type = "matrix")$D
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states + DD %*% dd
    } else {
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states
    }
    ## Var in model fits
    VtT <- MARSSkfss(MLEobj)$VtT
    VV <- NULL
    for (tt in 1:TT) {
      RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
                                                           , tt] %*% t(ZZ)
      SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
        t(MLEobj$states[, tt, drop = FALSE])
      VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
    }
    SE <- sqrt(VV)
    ## upper & lower (1-alpha)% CI
    fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
    fits$lo <- qnorm(alpha/2) * SE + fits$ex
    return(fits)
  }
  
  
  ## get model fits & CI's
  mod_fit <- get_DFA_fits(dfa_1)
  ## plot the fits
  ylbl <- obs_data
  par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
                                                               0, 0, 0))
  for (i in 1:N_ts) {
    up <- mod_fit$up[i, ]
    mn <- mod_fit$ex[i, ]
    lo <- mod_fit$lo[i, ]
    plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
         cex.lab = 1.2, ylim = c(min(lo), max(up)))
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
    lines(w_ts, dat[i, ], pch = 16, col = clr[i], lwd = 2)
    lines(w_ts, up, col = "darkgray")
    lines(w_ts, mn, col = "black", lwd = 2)
    lines(w_ts, lo, col = "darkgray")
    
  }
  
  ############## Regression work
  
  tempNEP <- array(unlist(NEP), dim = c(1, dim(NEP)))
  
  # Fit regression model
  sat.mod1 <- lm(tempNEP[1:365] ~ proc_rot[1,1:365]) # regression formula
  #sat.mod2 <- lm(tempNEP[1:365] ~ proc_rot[2,1:365]) # regression formula
  
  
  ################# SAVING DATA EACH YEAR
  Z_est_T4[,,z] <- Z_est
  lm_r2_T4[1,1,z] <- summary(sat.mod1)$r.squared
  lm_coef_T4[1,1,z] <- summary(sat.mod1)$coefficients[2, 1]
  lm_coef_T4[1,2,z] <- summary(sat.mod1)$coefficients[1, 1]
}
######## END OF T4 LOOP






######## START OF T5 LOOP
for (z in 1:8){
  year<-2008+z
  ## use only the 10 years from 1980-1989
  yr_frst <- year
  yr_last <- year
  mRUE_dat <- all_dat[all_dat[, "Year"] >= yr_frst & all_dat[, 
                                                              "Year"] <= yr_last, ]
  ## create vector of group names
  obs_data <- c("T5_day_WUEs", "T5_day_LUEs", "T5_day_CUEs")
  
  
  ## get only the RUEs
  dat_working <- mRUE_dat[, obs_data]
  NEP <- mRUE_dat[,"T5_day_NEP"]
  
  ## transpose data so time goes across columns
  dat_working <- t(dat_working)
  
  #### NAN counting
  nan_count_T5[1,1,z] = sum(is.na(dat_working[1,]))
  nan_count_T5[2,1,z] = sum(is.na(dat_working[2,]))
  nan_count_T5[3,1,z] = sum(is.na(dat_working[3,]))
  
  ## get number of time series
  N_ts <- dim(dat_working)[1]
  ## get length of time series
  TT <- dim(dat_working)[2]
  
  y_bar <- apply(dat_working, 1, mean, na.rm = TRUE)
  dat <- dat_working - y_bar
  rownames(dat) <- rownames(dat_working)
  
  ## number of processes
  mm <- 1
  ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
  BB <- "identity"  # diag(mm)
  ## 'uu' is a column vector of 0's
  uu <- "zero"  # matrix(0,mm,1)
  ## 'CC' and 'cc' are for covariates
  CC <- "zero"  # matrix(0,mm,1)
  cc <- "zero"  # matrix(0,1,wk_last)
  ## 'QQ' is identity
  QQ <- "identity"  # diag(mm)
  
  
  ## 'ZZ' is loadings matrix, 1x1
  Z_vals <- list("z11", "z21", "z31")
  ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 1, byrow = TRUE)
  ZZ
  
  ## 'aa' is the offset/scaling
  aa <- "zero"
  ## 'DD' and 'd' are for covariates
  DD <- "zero"  # matrix(0,mm,1)
  dd <- "zero"  # matrix(0,1,wk_last)
  ## 'RR' is var-cov matrix for obs errors
  RR <- "diagonal and unequal"
  
  ## list with specifications for model vectors/matrices
  mod_list <- list(B = BB, U = uu, C = CC, c = cc, Q = QQ, Z = ZZ, 
                   A = aa, D = DD, d = dd, R = RR)
  ## list with model inits
  init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
  ## list with model control parameters
  con_list <- list(maxit = 10000, allow.degen = TRUE)
  
  ## fit MARSS
  dfa_1 <- MARSS(y = dat, model = mod_list, inits = init_list, 
                 control = con_list)
  
  ## get the estimated ZZ
  Z_est <- coef(dfa_1, type = "matrix")$Z
  
  ## get the inverse of the rotation matrix
  #H_inv <- varimax(Z_est)$rotmat
  
  ## rotate factor loadings
  #Z_rot = Z_est %*% H_inv
  ## rotate processes
  proc_rot = Z_est %*% dfa_1$states
  
  
  
  ###### two stage loading plots
  ylbl <- obs_data
  w_ts <- seq(dim(dat)[2])
  layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
  ## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
  par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))
  ## plot the processes
  for (i in 1:mm) {
    ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
    ## set up plot area
    plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
         xlab = "", ylab = "", xaxt = "n")
    ## draw zero-line
    abline(h = 0, col = "gray")
    ## plot trend line
    lines(w_ts, proc_rot[i, ], lwd = 2)
    lines(w_ts, proc_rot[i, ], lwd = 2)
    ## add panel labels
    mtext(paste("State", i), side = 3, line = 0.5)
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
  }
  ## plot the loadings
  minZ <- 0
  ylm <- c(-1, 1) * max(abs(Z_est))
  for (i in 1:mm) {
    plot(c(1:N_ts)[abs(Z_est[, i]) > minZ], as.vector(Z_est[abs(Z_est[, 
                                                                      i]) > minZ, i]), type = "h", lwd = 2, xlab = "", ylab = "", 
         xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col = clr)
    for (j in 1:N_ts) {
      if (Z_est[j, i] > minZ) {
        text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
             col = clr[j])
      }
      if (Z_est[j, i] < -minZ) {
        text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
             col = clr[j])
      }
      abline(h = 0, lwd = 1.5, col = "gray")
    }
    mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
  }
  
  
  
  get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
    ## empty list for results
    fits <- list()
    ## extra stuff for var() calcs
    Ey <- MARSS:::MARSShatyt(MLEobj)
    ## model params
    ZZ <- coef(MLEobj, type = "matrix")$Z
    ## number of obs ts
    nn <- dim(Ey$ytT)[1]
    ## number of time steps
    TT <- dim(Ey$ytT)[2]
    ## get the inverse of the rotation matrix
    #H_inv <- varimax(ZZ)$rotmat
    ## check for covars
    if (!is.null(dd)) {
      DD <- coef(MLEobj, type = "matrix")$D
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states + DD %*% dd
    } else {
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states
    }
    ## Var in model fits
    VtT <- MARSSkfss(MLEobj)$VtT
    VV <- NULL
    for (tt in 1:TT) {
      RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
                                                           , tt] %*% t(ZZ)
      SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
        t(MLEobj$states[, tt, drop = FALSE])
      VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
    }
    SE <- sqrt(VV)
    ## upper & lower (1-alpha)% CI
    fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
    fits$lo <- qnorm(alpha/2) * SE + fits$ex
    return(fits)
  }
  
  
  ## get model fits & CI's
  mod_fit <- get_DFA_fits(dfa_1)
  ## plot the fits
  ylbl <- obs_data
  par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
                                                               0, 0, 0))
  for (i in 1:N_ts) {
    up <- mod_fit$up[i, ]
    mn <- mod_fit$ex[i, ]
    lo <- mod_fit$lo[i, ]
    plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
         cex.lab = 1.2, ylim = c(min(lo), max(up)))
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
    lines(w_ts, dat[i, ], pch = 16, col = clr[i], lwd = 2)
    lines(w_ts, up, col = "darkgray")
    lines(w_ts, mn, col = "black", lwd = 2)
    lines(w_ts, lo, col = "darkgray")
    
  }
  
  ############## Regression work
  
  tempNEP <- array(unlist(NEP), dim = c(1, dim(NEP)))
  
  # Fit regression model
  sat.mod1 <- lm(tempNEP[1:365] ~ proc_rot[1,1:365]) # regression formula
  #sat.mod2 <- lm(tempNEP[1:365] ~ proc_rot[2,1:365]) # regression formula
  
  
  ################# SAVING DATA EACH YEAR
  Z_est_T5[,,z] <- Z_est
  lm_r2_T5[1,1,z] <- summary(sat.mod1)$r.squared
  lm_coef_T5[1,1,z] <- summary(sat.mod1)$coefficients[2, 1]
  lm_coef_T5[1,2,z] <- summary(sat.mod1)$coefficients[1, 1]
}
######## END OF T5 LOOP






######## START OF T6 LOOP
for (z in 1:8){
  year<-2008+z
  ## use only the 10 years from 1980-1989
  yr_frst <- year
  yr_last <- year
  mRUE_dat <- all_dat[all_dat[, "Year"] >= yr_frst & all_dat[, 
                                                              "Year"] <= yr_last, ]
  ## create vector of group names
  obs_data <- c("T6_day_WUEs", "T6_day_LUEs", "T6_day_CUEs")
  
  
  ## get only the RUEs
  dat_working <- mRUE_dat[, obs_data]
  NEP <- mRUE_dat[,"T6_day_NEP"]
  
  ## transpose data so time goes across columns
  dat_working <- t(dat_working)
  
  #### NAN counting
  nan_count_T6[1,1,z] = sum(is.na(dat_working[1,]))
  nan_count_T6[2,1,z] = sum(is.na(dat_working[2,]))
  nan_count_T6[3,1,z] = sum(is.na(dat_working[3,]))
  
  ## get number of time series
  N_ts <- dim(dat_working)[1]
  ## get length of time series
  TT <- dim(dat_working)[2]
  
  y_bar <- apply(dat_working, 1, mean, na.rm = TRUE)
  dat <- dat_working - y_bar
  rownames(dat) <- rownames(dat_working)
  

  clr <- c("blue", "darkgreen", "darkred")

  
  ## number of processes
  mm <- 1
  ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
  BB <- "identity"  # diag(mm)
  ## 'uu' is a column vector of 0's
  uu <- "zero"  # matrix(0,mm,1)
  ## 'CC' and 'cc' are for covariates
  CC <- "zero"  # matrix(0,mm,1)
  cc <- "zero"  # matrix(0,1,wk_last)
  ## 'QQ' is identity
  QQ <- "identity"  # diag(mm)
  

  ## 'ZZ' is loadings matrix, 1x1
  Z_vals <- list("z11", "z21", "z31")
  ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 1, byrow = TRUE)
  ZZ
  
  ## 'aa' is the offset/scaling
  aa <- "zero"
  ## 'DD' and 'd' are for covariates
  DD <- "zero"  # matrix(0,mm,1)
  dd <- "zero"  # matrix(0,1,wk_last)
  ## 'RR' is var-cov matrix for obs errors
  RR <- "diagonal and unequal"
  
  ## list with specifications for model vectors/matrices
  mod_list <- list(B = BB, U = uu, C = CC, c = cc, Q = QQ, Z = ZZ, 
                   A = aa, D = DD, d = dd, R = RR)
  ## list with model inits
  init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
  ## list with model control parameters
  con_list <- list(maxit = 10000, allow.degen = TRUE)
  
  ## fit MARSS
  dfa_1 <- MARSS(y = dat, model = mod_list, inits = init_list, 
                 control = con_list)
  
  ## get the estimated ZZ
  Z_est <- coef(dfa_1, type = "matrix")$Z
  
  ## get the inverse of the rotation matrix
  #H_inv <- varimax(Z_est)$rotmat
  
  ## rotate factor loadings
  #Z_rot = Z_est %*% H_inv
  ## rotate processes
  proc_rot = Z_est %*% dfa_1$states
  
  
  
  ###### two stage loading plots
  ylbl <- obs_data
  w_ts <- seq(dim(dat)[2])
  layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
  ## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
  par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))
  ## plot the processes
  for (i in 1:mm) {
    ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
    ## set up plot area
    plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
         xlab = "", ylab = "", xaxt = "n")
    ## draw zero-line
    abline(h = 0, col = "gray")
    ## plot trend line
    lines(w_ts, proc_rot[i, ], lwd = 2)
    lines(w_ts, proc_rot[i, ], lwd = 2)
    ## add panel labels
    mtext(paste("State", i), side = 3, line = 0.5)
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
  }
  ## plot the loadings
  minZ <- 0
  ylm <- c(-1, 1) * max(abs(Z_est))
  for (i in 1:mm) {
    plot(c(1:N_ts)[abs(Z_est[, i]) > minZ], as.vector(Z_est[abs(Z_est[, 
                                                                      i]) > minZ, i]), type = "h", lwd = 2, xlab = "", ylab = "", 
         xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col = clr)
    for (j in 1:N_ts) {
      if (Z_est[j, i] > minZ) {
        text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
             col = clr[j])
      }
      if (Z_est[j, i] < -minZ) {
        text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
             col = clr[j])
      }
      abline(h = 0, lwd = 1.5, col = "gray")
    }
    mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
  }
  
  
  
  
  get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
    ## empty list for results
    fits <- list()
    ## extra stuff for var() calcs
    Ey <- MARSS:::MARSShatyt(MLEobj)
    ## model params
    ZZ <- coef(MLEobj, type = "matrix")$Z
    ## number of obs ts
    nn <- dim(Ey$ytT)[1]
    ## number of time steps
    TT <- dim(Ey$ytT)[2]
    ## get the inverse of the rotation matrix
    #H_inv <- varimax(ZZ)$rotmat
    ## check for covars
    if (!is.null(dd)) {
      DD <- coef(MLEobj, type = "matrix")$D
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states + DD %*% dd
    } else {
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states
    }
    ## Var in model fits
    VtT <- MARSSkfss(MLEobj)$VtT
    VV <- NULL
    for (tt in 1:TT) {
      RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
                                                           , tt] %*% t(ZZ)
      SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
        t(MLEobj$states[, tt, drop = FALSE])
      VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
    }
    SE <- sqrt(VV)
    ## upper & lower (1-alpha)% CI
    fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
    fits$lo <- qnorm(alpha/2) * SE + fits$ex
    return(fits)
  }
  
  
  ## get model fits & CI's
  mod_fit <- get_DFA_fits(dfa_1)
  ## plot the fits
  ylbl <- obs_data
  par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
                                                               0, 0, 0))
  for (i in 1:N_ts) {
    up <- mod_fit$up[i, ]
    mn <- mod_fit$ex[i, ]
    lo <- mod_fit$lo[i, ]
    plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
         cex.lab = 1.2, ylim = c(min(lo), max(up)))
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
    lines(w_ts, dat[i, ], pch = 16, col = clr[i], lwd = 2)
    lines(w_ts, up, col = "darkgray")
    lines(w_ts, mn, col = "black", lwd = 2)
    lines(w_ts, lo, col = "darkgray")
    
  }
  
  ############## Regression work
  
  tempNEP <- array(unlist(NEP), dim = c(1, dim(NEP)))
  
  # Fit regression model
  sat.mod1 <- lm(tempNEP[1:365] ~ proc_rot[1,1:365]) # regression formula
  #sat.mod2 <- lm(tempNEP[1:365] ~ proc_rot[2,1:365]) # regression formula
  
  
  ################# SAVING DATA EACH YEAR
  Z_est_T6[,,z] <- Z_est
  lm_r2_T6[1,1,z] <- summary(sat.mod1)$r.squared
  lm_coef_T6[1,1,z] <- summary(sat.mod1)$coefficients[2, 1]
  lm_coef_T6[1,2,z] <- summary(sat.mod1)$coefficients[1, 1]
}
######## END OF T6 LOOP






######## START OF T7 LOOP
for (z in 1:8){
  year<-2008+z
  ## use only the 10 years from 1980-1989
  yr_frst <- year
  yr_last <- year
  mRUE_dat <- all_dat[all_dat[, "Year"] >= yr_frst & all_dat[, 
                                                              "Year"] <= yr_last, ]
  ## create vector of group names
  obs_data <- c("T7_day_WUEs", "T7_day_LUEs", "T7_day_CUEs")
  
  
  ## get only the RUEs
  dat_working <- mRUE_dat[, obs_data]
  NEP <- mRUE_dat[,"T7_day_NEP"]
  
  ## transpose data so time goes across columns
  dat_working <- t(dat_working)
  
  #### NAN counting
  nan_count_T7[1,1,z] = sum(is.na(dat_working[1,]))
  nan_count_T7[2,1,z] = sum(is.na(dat_working[2,]))
  nan_count_T7[3,1,z] = sum(is.na(dat_working[3,]))
  
  ## get number of time series
  N_ts <- dim(dat_working)[1]
  ## get length of time series
  TT <- dim(dat_working)[2]
  
  y_bar <- apply(dat_working, 1, mean, na.rm = TRUE)
  dat <- dat_working - y_bar
  rownames(dat) <- rownames(dat_working)
  

  clr <- c("blue", "darkgreen", "darkred")

  
  ## number of processes
  mm <- 1
  ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
  BB <- "identity"  # diag(mm)
  ## 'uu' is a column vector of 0's
  uu <- "zero"  # matrix(0,mm,1)
  ## 'CC' and 'cc' are for covariates
  CC <- "zero"  # matrix(0,mm,1)
  cc <- "zero"  # matrix(0,1,wk_last)
  ## 'QQ' is identity
  QQ <- "identity"  # diag(mm)
  
  
  ## 'ZZ' is loadings matrix, 1x1
  Z_vals <- list("z11", "z21", "z31")
  ZZ <- matrix(Z_vals, nrow = N_ts, ncol = 1, byrow = TRUE)
  ZZ
  
  ## 'aa' is the offset/scaling
  aa <- "zero"
  ## 'DD' and 'd' are for covariates
  DD <- "zero"  # matrix(0,mm,1)
  dd <- "zero"  # matrix(0,1,wk_last)
  ## 'RR' is var-cov matrix for obs errors
  RR <- "diagonal and unequal"
  
  ## list with specifications for model vectors/matrices
  mod_list <- list(B = BB, U = uu, C = CC, c = cc, Q = QQ, Z = ZZ, 
                   A = aa, D = DD, d = dd, R = RR)
  ## list with model inits
  init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
  ## list with model control parameters
  con_list <- list(maxit = 10000, allow.degen = TRUE)
  
  ## fit MARSS
  dfa_1 <- MARSS(y = dat, model = mod_list, inits = init_list, 
                 control = con_list)
  
  ## get the estimated ZZ
  Z_est <- coef(dfa_1, type = "matrix")$Z
  
  ## get the inverse of the rotation matrix
  #H_inv <- varimax(Z_est)$rotmat
  
  ## rotate factor loadings
  #Z_rot = Z_est %*% H_inv
  ## rotate processes
  proc_rot = Z_est %*% dfa_1$states
  
  
  
  ###### two stage loading plots
  ylbl <- obs_data
  w_ts <- seq(dim(dat)[2])
  layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
  ## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
  par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))
  ## plot the processes
  for (i in 1:mm) {
    ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
    ## set up plot area
    plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
         xlab = "", ylab = "", xaxt = "n")
    ## draw zero-line
    abline(h = 0, col = "gray")
    ## plot trend line
    lines(w_ts, proc_rot[i, ], lwd = 2)
    lines(w_ts, proc_rot[i, ], lwd = 2)
    ## add panel labels
    mtext(paste("State", i), side = 3, line = 0.5)
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
  }
  ## plot the loadings
  minZ <- 0
  ylm <- c(-1, 1) * max(abs(Z_est))
  for (i in 1:mm) {
    plot(c(1:N_ts)[abs(Z_est[, i]) > minZ], as.vector(Z_est[abs(Z_est[, 
                                                                      i]) > minZ, i]), type = "h", lwd = 2, xlab = "", ylab = "", 
         xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col = clr)
    for (j in 1:N_ts) {
      if (Z_est[j, i] > minZ) {
        text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
             col = clr[j])
      }
      if (Z_est[j, i] < -minZ) {
        text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
             col = clr[j])
      }
      abline(h = 0, lwd = 1.5, col = "gray")
    }
    mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
  }
  
  
  
  
  get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
    ## empty list for results
    fits <- list()
    ## extra stuff for var() calcs
    Ey <- MARSS:::MARSShatyt(MLEobj)
    ## model params
    ZZ <- coef(MLEobj, type = "matrix")$Z
    ## number of obs ts
    nn <- dim(Ey$ytT)[1]
    ## number of time steps
    TT <- dim(Ey$ytT)[2]
    ## get the inverse of the rotation matrix
    #H_inv <- varimax(ZZ)$rotmat
    ## check for covars
    if (!is.null(dd)) {
      DD <- coef(MLEobj, type = "matrix")$D
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states + DD %*% dd
    } else {
      ## model expectation
      fits$ex <- ZZ %*% MLEobj$states
    }
    ## Var in model fits
    VtT <- MARSSkfss(MLEobj)$VtT
    VV <- NULL
    for (tt in 1:TT) {
      RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
                                                           , tt] %*% t(ZZ)
      SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
        t(MLEobj$states[, tt, drop = FALSE])
      VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
    }
    SE <- sqrt(VV)
    ## upper & lower (1-alpha)% CI
    fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
    fits$lo <- qnorm(alpha/2) * SE + fits$ex
    return(fits)
  }
  
  
  ## get model fits & CI's
  mod_fit <- get_DFA_fits(dfa_1)
  ## plot the fits
  ylbl <- obs_data
  par(mfrow = c(N_ts, 1), mai = c(0.5, 0.7, 0.1, 0.1), omi = c(0, 
                                                               0, 0, 0))
  for (i in 1:N_ts) {
    up <- mod_fit$up[i, ]
    mn <- mod_fit$ex[i, ]
    lo <- mod_fit$lo[i, ]
    plot(w_ts, mn, xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", 
         cex.lab = 1.2, ylim = c(min(lo), max(up)))
    axis(1, 365 * (0:dim(dat_working)[2]) + 1, yr_frst + 0:dim(dat_working)[2])
    lines(w_ts, dat[i, ], pch = 16, col = clr[i], lwd = 2)
    lines(w_ts, up, col = "darkgray")
    lines(w_ts, mn, col = "black", lwd = 2)
    lines(w_ts, lo, col = "darkgray")
    
  }
  
  ############## Regression work
  
  tempNEP <- array(unlist(NEP), dim = c(1, dim(NEP)))
  
  # Fit regression model
  sat.mod1 <- lm(tempNEP[1:365] ~ proc_rot[1,1:365]) # regression formula
  #sat.mod2 <- lm(tempNEP[1:365] ~ proc_rot[2,1:365]) # regression formula
  
  ################# SAVING DATA EACH YEAR
  Z_est_T7[,,z] <- Z_est
  lm_r2_T7[1,1,z] <- summary(sat.mod1)$r.squared
  lm_coef_T7[1,1,z] <- summary(sat.mod1)$coefficients[2, 1]
  lm_coef_T7[1,2,z] <- summary(sat.mod1)$coefficients[1, 1]
}
######## END OF T7 LOOP







all_nan_count = array(, dim=c(3,8,7))
all_nan_count[,,1] = nan_count_T1[,1,]
all_nan_count[,,2] = nan_count_T2[,1,]
all_nan_count[,,3] = nan_count_T3[,1,]
all_nan_count[,,4] = nan_count_T4[,1,]
all_nan_count[,,5] = nan_count_T5[,1,]
all_nan_count[,,6] = nan_count_T6[,1,]
all_nan_count[,,7] = nan_count_T7[,1,]



all_loadingfactors = array(, dim=c(3,8,7))

all_loadingfactors[,,1] = Z_est_T1[,1,]
all_loadingfactors[,,2] = Z_est_T2[,1,]
all_loadingfactors[,,3] = Z_est_T3[,1,]
all_loadingfactors[,,4] = Z_est_T4[,1,]  
all_loadingfactors[,,5] = Z_est_T5[,1,]  
all_loadingfactors[,,6] = Z_est_T6[,1,]  
all_loadingfactors[,,7] = Z_est_T7[,1,]  

mean_loadingfactors = apply(all_loadingfactors,c(1:2),mean)
sd_loadingfactors = apply(all_loadingfactors,c(1:2),sd)


all_lm_r2 = array(, dim=c(8,7))

all_lm_r2[,1] = lm_r2_T1[,1,]
all_lm_r2[,2] = lm_r2_T2[,1,]
all_lm_r2[,3] = lm_r2_T3[,1,]
all_lm_r2[,4] = lm_r2_T4[,1,]  
all_lm_r2[,5] = lm_r2_T5[,1,]  
all_lm_r2[,6] = lm_r2_T6[,1,]  
all_lm_r2[,7] = lm_r2_T7[,1,] 


all_lm_coef = array(, dim=c(2,8,7))

all_lm_coef[,,1] = lm_coef_T1[1,,]
all_lm_coef[,,2] = lm_coef_T2[1,,]
all_lm_coef[,,3] = lm_coef_T3[1,,]
all_lm_coef[,,4] = lm_coef_T4[1,,]  
all_lm_coef[,,5] = lm_coef_T5[1,,]  
all_lm_coef[,,6] = lm_coef_T6[1,,]  
all_lm_coef[,,7] = lm_coef_T7[1,,]  




plot(mean_loadingfactors[1,], type="b", xlab = "WUE", ylim=c(-0.2, 0.6), main="NEP stage")
plot(mean_loadingfactors[2,], type="b", xlab = "LUE", ylim=c(-0.2, 0.6))
plot(mean_loadingfactors[3,], type="b", xlab = "CUE", ylim=c(-0.2, 0.6))






xtime <- 2009:2016
#layout(matrix(c(1), mm, 2), widths = c(1))


plot.new()
plot.window(c(2008,2017),c(-0,0.7))
axis(1, at = 2008:2017)
axis(2)
lines(xtime,mean_loadingfactors[1,], col = "blue")
#segments(xtime,mean_loadingfactors[1,]-sd_loadingfactors[1,],xtime,mean_loadingfactors[1,]+sd_loadingfactors[1,], col = "blue")
text(2016.2,mean_loadingfactors[1,8],"WUE", col = "blue")

lines(xtime,mean_loadingfactors[2,], col = "darkgreen")
#segments(xtime,mean_loadingfactors[2,]-sd_loadingfactors[2,],xtime,mean_loadingfactors[2,]+sd_loadingfactors[2,], col = "darkgreen")
text(2016.2,mean_loadingfactors[2,8],"LUE", col = "darkgreen")

lines(xtime,mean_loadingfactors[3,], col = "darkred")
text(2016.2,mean_loadingfactors[3,8],"CUE", col = "darkred")

title(main="Stage 1")






#layout(matrix(c(2,1), mm, 2), widths = c(1))
plot.new()
plot.window(c(2008,2017),c(0,0.45))
axis(1, at = 2008:2017)
axis(2)
lines(xtime,mean_loadingfactors[1,], col = "blue")
text(2016.2,mean_loadingfactors[1,8],"WUE", col = "blue")
points(xtime,all_loadingfactors[1,,1])
#text(2016.2,all_loadingfactors[1,8,1],"T1")
points(xtime,all_loadingfactors[1,,2])
#text(2016.2,all_loadingfactors[1,8,2],"T2")
points(xtime,all_loadingfactors[1,,3])
#text(2016.2,all_loadingfactors[1,8,3],"T3")
points(xtime,all_loadingfactors[1,,4])
#text(2016.2,all_loadingfactors[1,8,4],"T4")
points(xtime,all_loadingfactors[1,,5])
#text(2016.2,all_loadingfactors[1,8,5],"T5")
points(xtime,all_loadingfactors[1,,6])
#text(2016.2,all_loadingfactors[1,8,6],"T6")
points(xtime,all_loadingfactors[1,,7])
#text(2016.2,all_loadingfactors[1,8,7],"T7")
title(main="WUE Loading Factors")





#layout(matrix(c(2,1), mm, 2), widths = c(1))
plot.new()
plot.window(c(2008,2017),c(0,1.2))
axis(1, at = 2008:2017)
axis(2)
lines(xtime,mean_loadingfactors[2,], col = "darkgreen")
text(2016.2,mean_loadingfactors[2,8],"LUE", col = "darkgreen")
points(xtime,all_loadingfactors[2,,1])
#text(2016.2,all_loadingfactors[2,8,1],"T1")
points(xtime,all_loadingfactors[2,,2])
#text(2016.2,all_loadingfactors[2,8,2],"T2")
points(xtime,all_loadingfactors[2,,3])
#text(2016.2,all_loadingfactors[2,8,3],"T3")
points(xtime,all_loadingfactors[2,,4])
#text(2016.2,all_loadingfactors[2,8,4],"T4")
points(xtime,all_loadingfactors[2,,5])
#text(2016.2,all_loadingfactors[2,8,5],"T5")
points(xtime,all_loadingfactors[2,,6])
#text(2016.2,all_loadingfactors[2,8,6],"T6")
points(xtime,all_loadingfactors[2,,7])
#text(2016.2,all_loadingfactors[2,8,7],"T7")
title(main="LUE Loading Factors")





#layout(matrix(c(2,1), mm, 2), widths = c(1))
plot.new()
plot.window(c(2008,2017),c(0,0.75))
axis(1, at = 2008:2017)
axis(2)
lines(xtime,mean_loadingfactors[3,], col = "darkred")
text(2016.2,mean_loadingfactors[3,8],"CUE", col = "darkred")
points(xtime,all_loadingfactors[3,,1])
#text(2016.2,all_loadingfactors[3,8,1],"T1")
points(xtime,all_loadingfactors[3,,2])
#text(2016.2,all_loadingfactors[3,8,2],"T2")
points(xtime,all_loadingfactors[3,,3])
#text(2016.2,all_loadingfactors[3,8,3],"T3")
points(xtime,all_loadingfactors[3,,4])
#text(2016.2,all_loadingfactors[3,8,4],"T4")
points(xtime,all_loadingfactors[3,,5])
#text(2016.2,all_loadingfactors[3,8,5],"T5")
points(xtime,all_loadingfactors[3,,6])
#text(2016.2,all_loadingfactors[3,8,6],"T6")
points(xtime,all_loadingfactors[3,,7])
#text(2016.2,all_loadingfactors[3,8,7],"T7")
title(main="CUE Loading Factors")






########## File output
write.csv(all_loadingfactors, file = "D:/desktop/a career is a lonelier path than a phd/KBS paper or papers/reprocessed data/all_loadingfactors_v2.dat")
write.csv(all_nan_count, file = "D:/desktop/a career is a lonelier path than a phd/KBS paper or papers/reprocessed data/all_nan_count_v2.dat")


###
write.csv(all_lm_r2, file = "D:/desktop/a career is a lonelier path than a phd/KBS paper or papers/reprocessed data/all_lm_r2_v2.dat")
write.csv(all_lm_coef, file = "D:/desktop/a career is a lonelier path than a phd/KBS paper or papers/reprocessed data/all_lm_coef_v2.dat")

#####
write.csv(tempNEP[1:365], file = "D:/desktop/a career is a lonelier path than a phd/KBS paper or papers/reprocessed data/regress_NEP_v2.dat")
write.csv(proc_rot[1,1:365], file = "D:/desktop/a career is a lonelier path than a phd/KBS paper or papers/reprocessed data/regress_model_v2.dat")

