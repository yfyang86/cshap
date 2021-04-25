quad_form <- function(A, x){
  # 二次型的快速计算方法
  return(sum(x * (A %*% x)))

}


ginverse=function(x, tol = sqrt(.Machine$double.eps))
{
  # return Moore-Penrose inverse of a matrix x
  if(length(dim(x)) > 2) stop("x must be a matrix or vector")
  svdX <- svd(x)
  if(is.complex(x))
    svdX$u <- Conj(svdX$u)
  NotZero <- svdX$d > tol * svdX$d[1]
  ans <- if(all(NotZero)) svdX$v %*% ((1/svdX$d) * t(svdX$
                                                       u)) else if(!any(NotZero)) {
                                                         if(is.matrix(x))
                                                           array(0, dim(x)[2:1])
                                                         else matrix(0, 1, length(x))
                                                       }
  else svdX$v[, NotZero] %*% ((1/svdX$d[NotZero]) * t(svdX$
                                                        u[, NotZero]))
  attr(ans, "rank") <- sum(NotZero)
  ans
}


AEM2r <- function(genodata, epsilon=1e-5, Nind = 50, p_initial = NULL, Allhaplo = NULL){
  # h p
  q = ncol(genodata)
  if (is.null(Allhaplo)) {Allhaplo <- haplo(q)}
  H = t(Allhaplo)
  r <- 2^q
  if (is.null(p_initial)||anyNA(p_initial)) {
    p = rep(1 / r, r)
  } else {
    p <- p_initial
  }
  # initial haplotype freq
  ii = 0
  truncthres <- 1 / r * 1e-3
  supp <- function(p){
    return(which(p > truncthres))
  }
  # Rh = rep(0, r)
  # Rh[supp(p)] <- 1
  # Rh <- Rh / length(supp(p))
  # IF <- rep(0, r)
  p = as.vector(p)
  repeat {
    ii = ii + 1
    #sigma0=var(genodata)
    #omega=apply(genodata,2,mean)
    # idx <- supp(p)
    idx <- which(p > 0)
    H2 <- H[idx,,drop = FALSE]
    p2 <- p[idx]
    omega <- colSums(H2 * p2)
    eta <- crossprod(H2, p2 * H2)
    sigma0 <- eta - tcrossprod(omega, omega)
    si <- ginverse(sigma0)
    IF <- rep(0, r)
    rh1res <- rep(0, length(idx))
    rh1expcoef <- rh1res
    rh2res <- matrix(0, nrow = length(idx), ncol = nrow(genodata))
    rh2expcoef <- rh2res
    for (j in seq_along(idx)) {
      h = H[idx[j], ]
      o_m_h <- omega - h
      rh1expcoef[j] <- -1 / (4 * Nind) * sum(o_m_h * (si %*% o_m_h))
      tmp = sweep(genodata / (2 * Nind), 2, omega, "-")

      tmpsi <- tmp %*% si
      expcoef <- - tmpsi %*% o_m_h - colSums(t(tmpsi) * t(tmp)) /2
      rh2expcoef[j,] <- expcoef
    }
    rh1expcoef <- rh1expcoef - max(rh1expcoef)
    rh2expcoef <- rh2expcoef - max(rh2expcoef)
    rh1res <- exp(rh1expcoef)
    rh2res <- exp(rh2expcoef)
    rh2meanres <- rowMeans(rh2res)
    rh2meanres <- rh2meanres/sum(rh2meanres)
    rh <- rh1res * rh2meanres
    rh <- rh / sum(rh)
    if (any(is.nan(rh))){
      return(NA)
    }
    IF[idx] <- rh
    Rh <- IF
    p.new = Rh * p
    p.new <- p.new / sum(p.new)
    p.new[p.new < truncthres] = 0
    p.new <- p.new / sum(p.new)
    delta <- sum(abs(p.new - p))
    # cat('i =', ii, ', delta =', delta, '\n')
    if (delta < epsilon | ii > 5000)
      break
    p <- p.new
  } #end of repeat
  idx <- which(p > 1e-5)
  H2 <- H[idx,,drop = FALSE]
  p2 <- p[idx]
  p2 <- p2 / sum(p2)
  return(list(freq = haplofreq0(q = q, ht = t(H2), p = p2)))
}

LAEM2r <- function(genodata, refHaplo = NULL, epsilon=1e-5, Nind = 50, p_initial = NULL){
  # h p
  q = ncol(genodata)
  if (is.null(refHaplo)) {refHaplo <- haplo(q)}
  H2 = t(refHaplo)
  r <- ncol(refHaplo)
  if (is.null(p_initial)||anyNA(p_initial)) {
    p = rep(1 / r, r)
  } else {
    p <- p_initial
  }
  # initial haplotype freq
  ii = 0
  truncthres <- 1 / r * 1e-3
  supp <- function(p){
    return(which(p > truncthres))
  }
  repeat {
    ii = ii + 1
    omega <- colSums(H2 * p)
    # eta <- crossprod(H2, diag(p2)) %*% H2
    eta <- crossprod(H2, p * H2)
    sigma0 <- eta - tcrossprod(omega, omega)
    si <- ginverse(sigma0)
    length(p) -> r
    rh1res <- rep(0, r)
    rh1expcoef <- rh1res
    rh2res <- matrix(0, nrow = r, ncol = nrow(genodata))
    rh2expcoef <- rh2res
    for (j in 1:r) {
      h = H2[j, ]
      o_m_h <- omega - h
      rh1expcoef[j] <- -1 / (4 * Nind) * sum(o_m_h * (si %*% o_m_h))
      tmp = sweep(genodata / (2 * Nind), 2, omega, "-")

      tmpsi <- tmp %*% si
      expcoef <- - tmpsi %*% o_m_h - colSums(t(tmpsi) * t(tmp)) /2
      rh2expcoef[j,] <- expcoef
    }
    rh1expcoef <- rh1expcoef - max(rh1expcoef)
    rh2expcoef <- rh2expcoef - max(rh2expcoef)
    rh1res <- exp(rh1expcoef)
    rh2res <- exp(rh2expcoef)
    rh2meanres <- rowMeans(rh2res)
    rh2meanres <- rh2meanres/sum(rh2meanres)
    rh <- rh1res * rh2meanres
    rh <- rh / sum(rh)
    if (any(is.nan(rh))){
      return(NA)
    }
    Rh <- rh
    p.new = Rh * p
    p.new <- p.new / sum(p.new)
    p.new[p.new < truncthres] = 0
    p.new <- p.new / sum(p.new)
    delta <- sum(abs(p.new - p))
    # cat('i =', ii, ', delta =', delta, '\n')
    if (delta < epsilon | ii > 5000)
      break
    p <- p.new
    S_idx <- p > 1e-5
    p <- p[S_idx]
    p <- p / sum(p)
    H2 <- H2[S_idx,,drop = FALSE]

  } #end of repeat

  return(list(freq = haplofreq0(q = q, ht = t(H2))))
} # end of AEM

AES2r <- function(genodata, rho, Nind = 50, p_initial = NULL, epsilon = 1e-6, Allhaplo = NULL) {
  if (abs(rho) < 1e-2) {
    return(
      AEM2r(genodata = genodata, epsilon = epsilon, Nind = Nind, p_initial = p_initial, Allhaplo = Allhaplo)
    )
  }
  q = ncol(genodata)
  T = nrow(genodata)
  if (is.null(Allhaplo)) {Allhaplo <- haplo(q)}
  H = t(Allhaplo)
  r <- 2^q
  if (is.null(p_initial)||anyNA(p_initial)) {
    p = rep(1 / r, r)
  } else {
    p <- p_initial }
  d = Inf
  truncthres <- min(0.1/r, max(1e-5, 1 / r * 1e-3))
  supp <- function(p){
    return(which(p > truncthres))
  }
  Rh <- rep(0, r)
  Rh[supp(p)] <- 1
  Rh <- Rh / length(supp(p))

  ii = 0
  repeat {
    ii = ii + 1
    idx <- supp(p)
    H2 <- H[idx,,drop = FALSE]
    p2 <- p[idx]
    omega <- colSums(H2 * p2)
    eta <- crossprod(H2, diag(p2)) %*% H2
    sigma0 <- eta - tcrossprod(omega, omega)
    sigma0 <- (1 + rho) * sigma0

    si = ginverse(sigma0)
    si1 = si / (2 * Nind - 2)
    si2 = si / (2 * Nind)

    IF <- rep(0, r) # importance factor
    for (j_idx in seq_along(idx)) {
      j <- idx[j_idx]
      h1 = H2[j_idx, ]
      Rj <- rep(0, T)
      for (i in 1:T)
      {
        Ai <- genodata[i, ]
        rw <- 0
        for (jj_idx in seq_along(idx))
        {
          jj <- idx[jj_idx]
          h2 = H2[jj_idx, ]
          pr <- (1 - rho) * p[jj] + rho * (j == jj)
          arg <- Ai - h1 - h2
          arg2 <- arg - 2 * (Nind - 1) * omega
          ep1 <- sum(arg2 * (si1 %*% arg2))
          arg3 <- Ai - 2 * Nind * omega
          ep2 <- sum(arg3 * (si2 %*% arg3))
          tjj <- exp(-(ep1 - ep2)/2)
          rw <- rw + pr * tjj

        }
        Rj[i] <- rw
      }
      IF[j] <- mean(Rj)
    } # end of all haplotypes
    Rh = IF
    Rh <- Rh / sum(Rh)
    p.new = Rh * p
    p.new[p.new < 0 ] <- 0
    tmp = sum(p.new)
    if (is.na(tmp)) return(NA)
    # in retrospect, this line should appear after the next line, but never mind
    p.new[p.new < 1 / 2 ^ q * 1e-3] = 0
    p.new = p.new / tmp
    delta <- sum(abs(p.new - p))
    # cat('i=', ii, ', delta=', delta, '\n')
    if (delta < epsilon | ii > 2000)
      break
    p = p.new
  } #end of repeat
  idx <- which(p > 1e-5)
  H2 <- H[idx,,drop = FALSE]
  p2 <- p[idx]
  p2 <- p2 / sum(p2)
  return(list(freq = haplofreq0(q = q, ht = t(H2), p = p2)))
} # end of AES

CSHAPAEM <- function(genodata, Nind = 50,
                     rho = 0, HWE = TRUE, epsilon=1e-5, .RM_Mat = NULL){
  q <- ncol(genodata)
  Allhaplo <- haplo(q)

  if (!HWE && !rho){
    ## Adaptive
    if (boot_rho(genodata, Nind = Nind, nboot = 100, check = FALSE)){
      rho <- estimate.rho(genodata, Nind = Nind, check = FALSE)
    } else {
      rho <- 0
    }
  }

  if (abs(rho) <= 1e-2){
    p_init <- CSHAP(genodata = genodata, Nind = Nind, rho = 0,quantile_thres = 0.5, HWE = TRUE,
                    Allhaplo = Allhaplo, .RM_Mat = .RM_Mat)
    res <- AEM2r(genodata, epsilon = epsilon, Nind = Nind, p_initial = p_init,
                 Allhaplo = Allhaplo)
  } else {
    p_init <- CSHAP(genodata = genodata, Nind = Nind, rho = rho, quantile_thres = 0.5,
                    Allhaplo = Allhaplo, .RM_Mat = .RM_Mat)
    res <- AES2r(genodata, rho = rho, Nind = Nind, p_initial = p_init,
                 epsilon = epsilon, Allhaplo = Allhaplo)
  }
  # HWD
  return(res)
}

ligation2LCSAEM <- function(genodata, ResA, ResB, rho = 0, Nind = 1, epsilon =1e-5){
  my_Haploset <- ligation_mat(ResA$freq$ht, ResB$freq$ht)
  q <- ncol(genodata)
  ref_p <- NA
  if (q <= 30 ) {
    ref_p <- LCSHAP0(genodata = genodata, rho = rho, my_Haploset = my_Haploset, Nind = Nind)
    ref_p <- as.vector(ref_p)
  }
  if (anyNA(ref_p)) {
    p <- outer(ResA$freq$p, ResB$freq$p)
    ref_p <- c(t(p))
  }
  S_idx <- ref_p > 0
  refHaplo <- my_Haploset[, S_idx, drop = FALSE]
  ref_p <- ref_p[S_idx]
  ref_p <- as.vector(ref_p)
  ref_p <- ref_p / sum(ref_p)

  res <- LAEM2r(genodata = genodata, refHaplo = refHaplo, Nind = Nind, epsilon = epsilon, p_initial = ref_p)
  return(res)
}

PLCSHAPLAEM <- function(genodata, Nind = 50, rho = 0, HWE = TRUE, epsilon =  1e-5, max_block_size = 13){
  ## block_fun 指的是每一块上使用的函数，需要输出为haplofreq对象
  ## 于每一块上使用来fastEM（目前是有权的）
  ## 连接优先使用CS
  q <- ncol(genodata)
  if (q <= 19) return(
    CSHAPAEM(genodata = genodata, rho = rho, HWE = HWE, Nind = Nind, epsilon = epsilon)
  )
  Pools <- nrow(genodata)
  pconfig <- partitionconfig(q, max_block_size = max_block_size)
  blocks <- nrow(pconfig)
  genodata_list <- vector(mode = 'list', length = blocks)
  names(genodata_list) <- paste0('block', 1:blocks)
  P_Res <-  genodata_list

  # if (! memoise::is.memoised(haplo)) haplo <-memoise::memoise(haplo)

  Allhaploi <-  haplo(pconfig[1,3])
  .RM_Mati <- haplo2RM(Allhaploi)

  for (i in seq_len(blocks)){
    genodata_list[[i]] <- genodata[, pconfig[i,1] : pconfig[i,2], drop = FALSE]
    if (i <= blocks - 2) {
      resi <- CSHAPAEM(genodata = genodata_list[[i]], rho = rho, HWE = HWE,
                        Nind = Nind, epsilon = epsilon, .RM_Mat = .RM_Mati)
    } else {
      resi <- CSHAPAEM(genodata = genodata_list[[i]], rho = rho, HWE = HWE,
                        Nind = Nind, epsilon = epsilon, .RM_Mat = NULL)
    }
    P_Res[[i]] <- resi
  }

  # return(pconfig)
  # names(P_Res) <- 1:blocks
  ### 对每一块进行计算
  while(blocks > 1){
    # L_step <- LigationStep(blocks)
    newblocks <- ceiling(blocks / 2)
    odds <- blocks %% 2
    new_p_len <- blocks %/% 2
    newpconfig <- matrix(0, ncol = 3, nrow = newblocks)
    colnames(newpconfig) <- colnames(pconfig)
    New_P_Res <- vector('list', newblocks)
    names(New_P_Res) <- 1:newblocks
    New_genodata_list <- New_ugdata_list <- New_P_Res
    for (i in seq_len(new_p_len)){
      ## 连接 newpconfig
      newpconfig[i, 1] <- newstart <- pconfig[2*i-1, 1]
      newpconfig[i, 2] <- newend <- pconfig[2*i, 2]
      newpconfig[i, 3] <- pconfig[2*i-1, 3]  + pconfig[2*i, 3]
      new_genodata <- genodata[, newstart : newend, drop = FALSE]
      # stopifnot( !l2dist(new_ugdata$udata[new_ugdata$uidx,], new_genodata) ) # 调试用，待会删
      new_out <- ligation2LCSAEM(genodata = new_genodata,  ResA = P_Res[[ 2 * i - 1]],  ResB = P_Res[[ 2 * i ]], rho = 0,
                                 Nind = Nind, epsilon = epsilon)

      New_P_Res[[i]] <- new_out
      New_genodata_list[[i]] <- new_genodata
    }
    if (odds) {
      New_P_Res[[new_p_len + 1]] <- last(P_Res)
      New_genodata_list[[new_p_len + 1]]<- last(genodata_list)
      newpconfig[newblocks, ] <- pconfig[blocks,]
    }
    blocks <- newblocks
    pconfig <- newpconfig
    P_Res <- New_P_Res
    genodata_list <- New_genodata_list
  }

  out <- New_P_Res[[1]]
  return(out)
}

