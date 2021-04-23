# Base functions ----------------------------------------------------------
# X1的第一列(的向量)和X2连接
ligation1 <- function(vec1, mat2){
  vec1 <- c(vec1)
  q2 <- nrow(mat2)
  n2 <- ncol(mat2)
  q1 <- length(vec1)
  A <- matrix(0, nrow = q1 + q2, ncol = n2)
  A[1:q1,] <- vec1
  A[(q1 + 1) : (q1 + q2),]  <- mat2
  return(A)
}
# 一个矩阵与矩阵连接
ligation_mat <- function(mat1, mat2){
  tmpres <- lapply(as.data.frame(mat1), ligation1, mat2 = mat2)
  return(do.call(cbind, tmpres))
}

# 两个频率相乘
ligation_p <- function(p1, p2){
  p <- outer(p1, p2)
  return(c(t(p)))
}


ligation.haploset0 <- function(X1, X2){
  ht <- ligation_mat(X1$ht, X2$ht)
  # htc <- outer(X1$htc, X2$htc, FUN = str_c)
  q <- X1$q + X2$q
  # htc <- c(t(htc))
  res <- list(q = q, ht = ht)
  class(res) <-'haploset0'
  return(res)
}

ligation.haploset <- function(X1, X2){
  ht <- ligation_mat(X1$ht, X2$ht)
  htc <- outer(X1$htc, X2$htc, FUN = str_c)
  q <- X1$q + X2$q
  htc <- c(t(htc))
  res <- list(q = q, ht = ht, htc = htc)
  class(res) <-'haploset'
  return(res)
}

ligation.haplofreq0 <- function(X1, X2){
  res <- ligation.haploset0(X1, X2)
  p <- outer(X1$p, X2$p)
  res$p <- c(t(p))
  class(res) <- c('haplofreq0','haploset0')
  return(res)
}

ligation.haplofreq <- function(X1, X2){
  res <- ligation.haploset(X1, X2)
  p <- outer(X1$p, X2$p)
  res$p <- c(t(p))
  class(res) <- c('haplofreq','haploset')
  return(res)
}

ligation <- function(x,y) UseMethod("ligation")


HtPhase_E_Step <- function(PhaseA, PhaseB, S_B, ref_p = NULL){
  ## 第一步，先从分布Phase得到连接后的ht的E_Step，由于是第一步，这里先不截断
  ## 之后可以当做 HF EM 来做
  nA <- nrow(PhaseA)
  nB <- nrow(PhaseB)
  
  homoA <- PhaseA[1, 1] == PhaseA[1, 2]
  homoB <- PhaseB[1, 1] == PhaseB[1, 2]
  
  if (homoA && homoB){
    ## A,B 均为纯合子
    Prob <- 1
    res <- cbind((PhaseA[1, 1] - 1) * S_B + PhaseB[1, 1:2, drop = FALSE], Prob)
    return(res)
  }
  
  if (homoA && !homoB){
    ## A 纯合子，B杂合子
    res <- (PhaseA[1, 1] - 1) * S_B + PhaseB[, 1:2, drop = FALSE]
    
    if (!is.null(ref_p)) {
      Prob <- ref_p[res[,2]] / sum(ref_p[res[,2]])
      if (anyNA(Prob)) Prob <- PhaseB[, 3] / sum(PhaseB[, 3])
    } else {
      Prob <- PhaseB[, 3] / sum(PhaseB[, 3])
    }
    
    res <- cbind(res, Prob)
    return(res[Prob > 0, , drop = FALSE])
  }
  
  if (!homoA && homoB){
    ## A为杂合子，B 为纯合子
    Prob <- PhaseA[, 3]
    Prob <- Prob / sum(Prob)
    res <- (PhaseA[, 1:2, drop = FALSE] - 1) * S_B + PhaseB[1, 1]
    
    if (!is.null(ref_p)) {
      Prob <- ref_p[res[,1]] / sum(ref_p[res[,1]])
      if (anyNA(Prob)) Prob <- PhaseA[, 3] / sum(PhaseA[, 3])
    } else {
      Prob <- PhaseA[, 3] / sum(PhaseA[, 3])
    }
    
    res <- cbind(res, Prob)
    return(res[Prob > 0, , drop = FALSE])
  }
  
  if (nA * nB == 1){
    ## A,B的相型均已经确定
    res <-  cbind(hAB1 = (PhaseA[1, 1] - 1) * S_B + PhaseB[1:2],
                  hAB2 = (PhaseA[1, 2] - 1) * S_B + PhaseB[2:1])
    
    if (!is.null(ref_p)){
      Prob <- ref_p[res[,1]] * ref_p[res[,2]]
      min(Prob) / 10 -> epsi
      compati1 <- ref_p[res[,1]] > 0 
      compati2 <- ref_p[res[,2]] > 0 
      Prob[compati1 & !compati2] <- epsi
      Prob[compati2 & !compati1] <- epsi
      Prob <- Prob / sum(Prob)
      if (anyNA(Prob)) Prob <- c(0.5, 0.5)
    } else {
      Prob <- c(0.5, 0.5)
    }
    res <- cbind(res, Prob)
    return(res[Prob > 0, , drop = FALSE])
  }
  
  jointHtP <- matrix(0, ncol = 5, nrow = 2 * nA * nB, dimnames = list(c(), c('hA1', 'hB1', 'hA2', 'hB2', 'Prob')))
  ## hA1 ~ hB1 = h1; hA2 ~ hB2 = h2; ~ 表示捆绑连接
  Prob <- numeric(2 * nA * nB)
  for (i in 1:nA){
    jointHtP[1:(2 * nB) + 2 * (i - 1) * nB, 1] <- PhaseA[i, 1]
    jointHtP[1:(2 * nB) + 2 * (i - 1) * nB, 3] <- PhaseA[i, 2]
    Prob[1:(2 * nB) + 2 * (i - 1) * nB] <- rep(PhaseA[i, 3] * PhaseB[, 3] / 2, each = 2) ## 一次性写完概率
    for (j in 1:nB){
      crow <- 2 * (i - 1) * nB + 2 * j - 1# idx of row
      jointHtP[crow, c(2, 4)] <- PhaseB[j, 1:2]
      jointHtP[crow + 1, c(2, 4)] <- PhaseB[j, 2:1]
      # jointHtP[crow : (crow + 1), 5] <- PhaseA[i, 3] * PhaseB[j, 3] / 2
    }
  }
  ## Mat 转换为 id 的方法：
  ## Ai Bj = Mat[i,j] = (i - 1) * S_B + j
  
  res <- cbind((jointHtP[, c(1,3)] - 1) * S_B + jointHtP[, c(2,4)])
  
  if (!is.null(ref_p)){
    Prob2 <- ref_p[res[,1]] * ref_p[res[,2]]
    min(Prob2) / 10 -> epsi
    compati1 <- ref_p[res[,1]] > 0 
    compati2 <- ref_p[res[,2]] > 0 
    Prob2[compati1 & !compati2] <- epsi
    Prob2[compati2 & !compati1] <- epsi
    Prob2 <- Prob2 / sum(Prob2)
    if (anyNA(Prob2)) Prob2 <- Prob / sum(Prob)
  } else {
    Prob2 <- Prob / sum(Prob)
  }
  
  
  res <- cbind(res, Prob2)
  colnames(res) <- c('hAB1', 'hAB2', 'Prob')
  return(res[Prob2 > 0, , drop = FALSE])
}

HtPhase_E_Steps <- function(PhaseAs, PhaseBs, S_B, ref_p = NULL) {
  # Res <- mapply(HtPhase_E_Step, PhaseAs, PhaseBs, S_B = S_B, ref_p = ref_p, SIMPLIFY = FALSE)  
  Res <- map2(.x = PhaseAs, .y = PhaseBs, .f = HtPhase_E_Step, S_B = S_B, ref_p = ref_p)
  return(Res)
}

if (FALSE){
  microbenchmark(
    for (i in 1:length(PhaseAs)){
      # cat(i, '\n')
      HtPhase_E_Step(PhaseA =  PhaseAs[[i]],  PhaseB =  PhaseBs[[i]], S_B, ref_p = ref_p) -> Res[[i]]
    },
    map2(PhaseAs, PhaseBs, HtPhase_E_Step, S_B = S_B, ref_p = ref_p)
  )
}


ligation.udata <- function(U1, U2) {
  uidx3 <- cbind(U1$uidx, U2$uidx)
  uniquedata(uidx3, MARGIN = 1, detail = TRUE) -> uudata3
  ugdata3 <- cbind(U1$udata[uudata3$udata[,1],],
                   U2$udata[uudata3$udata[,2],])
  rownames(ugdata3) <- seq_len(nrow(ugdata3))
  colnames(ugdata3) <- c(colnames(U1$udata), colnames(U2$udata))
  udata3 <- uudata3
  udata3$udata <- ugdata3
  return(udata3)
}


fastligationEM  <- function(ResA, ResB, udataAB, ref_p = NULL, truncthres2 = 1e-3, uniqueDiplo_in = FALSE, uniqueDiplo_out = FALSE, transMat = FALSE) {
  ## 区别在于加上了udata，速度未必有提升，但是对于加权EM求lone.htn是非常重要的。
  ## uniqueDiplo_in指的是输入唯一解，uniqueDiplo_out指的是输出唯一解
  S_A <- length(ResA$freq$p)
  S_B <- length(ResB$freq$p)
  N_Inds <- length(ResA$HtPhase$Phase)
  out_hf <- Lig_hf <- ligation.haplofreq0(X1 = ResA$freq, X2 = ResB$freq)
  r <- S_A * S_B
  if (is.null(ref_p)) {
    freq <- ref_p <- Lig_hf$p 
  } else {
    freq <- ref_p  
  }
  matchidx <- sapply(udataAB$matchidx, function(o){o[[1]]})
  ucount <- udataAB$count
  
  if (!uniqueDiplo_in | is.null(is.null(ResA$HtDiplo))) {
    UE_Results <- HtPhase_E_Steps(PhaseAs = ResA$HtPhase$Phase[matchidx], PhaseBs = ResB$HtPhase$Phase[matchidx],
                                  S_B = S_B, ref_p = ref_p)
    
  } else {
    UE_Results <- HtDiplo_E_Steps(HtDAs = ResA$HtDiplo$Diplo[matchidx, , drop = FALSE], 
                                  HtDBs = ResB$HtDiplo$Diplo[matchidx, , drop = FALSE], S_B = S_B)
  }
  
  UM_Res <- M_Steps(E_Results = UE_Results, r = r)
  newfreq <- c(UM_Res %*% ucount) / N_Inds
  delta <- sum(abs(newfreq - freq))
  
  interestidx <- which(sapply(UE_Results, nrow) > 1)
  n_interest <- length(interestidx)
  UM_Res_fixed <- UM_Res[, -interestidx, drop = FALSE] %*% ucount[-interestidx]
  
  while (delta > 1e-6 & n_interest){
    # cat(delta, '\t', n_interest, '\n')
    freq <- newfreq
    UE_Results_interest <- UE_Results[interestidx]
    UE_Results_interest <- E_Steps_withdiplo(E_Results = UE_Results_interest, freq = freq, truncthres2 = truncthres2)
    UE_Results[interestidx] <- UE_Results_interest
    new_int_id <- sapply(UE_Results_interest, nrow) > 1 
    new_n_interest <- sum(new_int_id)
    UM_Res_interest <- M_Steps(UE_Results_interest, r = r)
    if (new_n_interest < n_interest) {
      UM_Res_fixed <- UM_Res_fixed + 
        UM_Res_interest[,!new_int_id, drop = FALSE] %*% ucount[interestidx][!new_int_id]
    }
    newfreq <- (UM_Res_interest[,new_int_id, drop = FALSE] %*% ucount[interestidx][new_int_id] + 
                  UM_Res_fixed) / N_Inds
    delta <- sum(abs(newfreq - freq))
    interestidx <- interestidx[new_int_id] 
    n_interest <- new_n_interest
  }
  
  newfreq <- c(newfreq / sum(newfreq))
  S_idx <- which(newfreq > 0)
  out_hf$p <- newfreq[S_idx]
  out_hf$ht <- Lig_hf$ht[, S_idx, drop = FALSE]
  out_hf$htc <- Lig_hf$htc[S_idx]
  
  Q_idx <- integer(S_A * S_B) ## 方便转换HtPhase到简化后的，这样可以不需要match
  Q_idx[S_idx] <- seq_along(S_idx)
  
  LE_Results <- UE_Results[udataAB$uidx]
  names(LE_Results) <- 1:N_Inds
  LHtPhase_list <- llply(.data = LE_Results, .fun = function(o){
    o[, 1] <- Q_idx[o[, 1]]
    o[, 2] <- Q_idx[o[, 2]]
    return(o)
  })  
  
  LHtPhase <- HtPhase(ht = out_hf, Phase = LHtPhase_list)
  # LHtDiplo <- as.HtDiplo.HtPhase(LHtPhase)
  result <- list(freq = out_hf, HtPhase = LHtPhase)
  
  if (transMat){
    transMat <- matrix(Q_idx, nrow = S_A, ncol = S_B, byrow = TRUE, dimnames = list(
      str_c('A', 1:S_A) , str_c('B', 1:S_B)
    )) 
    ### 这里暂时先用Matrix方便调试，待会需要提高性能的时候再换
    transMatP <- matrix(newfreq, nrow = S_A, ncol = S_B, byrow = TRUE, dimnames = list(
      str_c('A', 1:S_A) , str_c('B', 1:S_B)
    ))
    
    transMatcP <- transMatP/ResA$freq$p # cP 表示的是条件概率，P(hB|hA)，也是 hA -> hB 的马尔科夫转移概率
    
    result$transMat <- transMat
    result$transMatjP <- transMatP ## JointP 
    result$transMatcP <- transMatcP ## transP
  }
  
  if (uniqueDiplo_out){
    result$HtDiplo <- as.HtDiplo.HtPhase(X_HtP = result$HtPhase)
  }
  return(result)
}


LCSHAP0 <- function(genodata, Nind = 1, rho = 0, my_Haploset = NULL, udata = NULL){
  Pools <- nrow(genodata)
  q <- ncol(genodata)
  if(is.null(udata) && Nind == 1) udata <- uniquegenodata(genodata, MARGIN = 1) 
  
  if (is.null(my_Haploset)) {
    E_Res0 <-  apply(udata$udata, 1, E_Step0_htn)
    res0 <- unique(unlist(sapply(E_Res0, c)))
    my_Haploset <- as_haplo2(q = q, htn = res0)
  }
  
  tmpe_e_MAF <- estimate.MAF(genodata, Nind = Nind, rho = rho)
  e_b <- estimate.b(genodata, e.MAF = tmpe_e_MAF, Nind = Nind)
  tmpRM <- haplo2RM(my_Haploset)
  restmp0 <- spams.lasso(X = e_b, D = tmpRM$phi, return_reg_path = FALSE,
                         lambda1 = 0 , mode = 'L2ERROR', pos = TRUE)
  tmpresult <- restmp0/tmpRM$tau
  # idx <- tmpresult@i + 1
  # ht = my_Haploset[,idx, drop = FALSE]
  # tmpresult <- tmpresult[idx]
  p <- tmpresult/sum(tmpresult)
  return(p)
  # if (is.null(htmat)) return(freq2haplofreq(p)) else 
  # return(as.haplofreq(htmat, p))
}


partitionconfig <- function(n, max_block_size = 10) {
  if ( n <=  2 * max_block_size){
    block_re1 <- floor(n/2)
    block_re2 <-ceiling(n/2)
    block_length <- c(block_re1, block_re2)
    block_start <- c(0, block_re1) + 1
    block_end <- c(block_re1, n)
  } else {
    blocks <- n %/% max_block_size
    if (n %% max_block_size){
      fixed <-  max_block_size * (blocks - 1)
      block_re1 <- floor((n - fixed)/2)
      block_re2 <- ceiling((n - fixed)/2)
      
      block_length <-  c(rep(max_block_size, blocks - 1),block_re1 ,block_re2)
      block_start <- c((0:(blocks - 1)) * max_block_size, fixed + block_re1) + 1
      block_end <- c((1:(blocks - 1)) * max_block_size, fixed + block_re1, fixed + block_re1 +block_re2 )
    } else { ## 此时整除
      block_length <- max_block_size
      block_start <-  (1:blocks - 1) * max_block_size + 1
      block_end <- 1:blocks * max_block_size
    }
    
  }
  
  config  <- cbind(start = block_start, end = block_end, length = block_length)
  return(config)
}

ligation2LCSEM <- function(genodata, udataAB = NULL, ResA, ResB, rho = 0, truncthres2 = 1e-04){
  my_Haploset <- ligation_mat(ResA$freq$ht, ResB$freq$ht)
  
  if (is.null(udataAB)) udataAB <- uniquegenodata(genodata, detail = FALSE)
  ref_p <- LCSHAP0(genodata = genodata, rho = rho, my_Haploset = my_Haploset, udata = udataAB)
  ref_p <- as.vector(ref_p)
  if (anyNA(ref_p)) ref_p <- NULL
  # 
  res <- fastligationEM(ResA, ResB, udataAB = udataAB, ref_p = ref_p, truncthres2 = truncthres2,
                        uniqueDiplo_in = FALSE, uniqueDiplo_out = FALSE, transMat = FALSE)
  
  return(res)
}


PLCSHAPLEM <- function(genodata, rho = 0, HWE = TRUE, truncthres2 = 1e-3, max_block_size = 13, inferphase = FALSE){
  ## block_fun 指的是每一块上使用的函数，需要输出为haplofreq对象
  ## 于每一块上使用来fastEM（目前是有权的）
  ## 连接优先使用CS
  q <- ncol(genodata)
  i_missing <- any(genodata == 3) || anyNA(genodata)
  if (q <= 19) return(CSHAPEM(genodata = genodata, rho = rho, HWE = HWE, truncthres2 = truncthres2, inferphase = inferphase))
  N_Inds <- nrow(genodata)
  pconfig <- partitionconfig(q, max_block_size = max_block_size)
  blocks <- nrow(pconfig)
  genodata_list <- vector(mode = 'list', length = blocks)
  names(genodata_list) <- paste0('block', 1:blocks)
  P_Res <- ugdata_list <- genodata_list
  
  # if (! memoise::is.memoised(haplo)) haplo <-memoise::memoise(haplo)
  
  Allhaploi <-  haplo(pconfig[1,3])
  .RM_Mati <- haplo2RM(Allhaploi)
  
  for (i in seq_len(blocks)){
    genodata_list[[i]] <- genodata[, pconfig[i,1] : pconfig[i,2], drop = FALSE]
    ugdata_list[[i]] <- uniquegenodata(genodata = genodata_list[[i]], MARGIN = 1, detail = TRUE, missing = i_missing)
    if (i <= blocks - 2) {
      resi <- CSHAPEM(genodata = genodata_list[[i]], rho = rho, HWE = HWE,
                      Allhaplo = Allhaploi, truncthres2 = truncthres2,
                         udata = ugdata_list[[i]], inferphase = TRUE, .RM_Mat = .RM_Mati)
    } else {
      resi <- CSHAPEM(genodata = genodata_list[[i]], rho = rho, HWE = HWE, 
                      Allhaplo = NULL, truncthres2 = truncthres2,
                         udata = ugdata_list[[i]], inferphase = TRUE, .RM_Mat = NULL)
    }
    P_Res[[i]] <- resi
  }
  
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
      
      new_ugdata <- ligation.udata(U1 = ugdata_list[[2 * i - 1]],
                                   U2 = ugdata_list[[2 * i]])
      new_genodata <- genodata[, newstart : newend, drop = FALSE]
      # stopifnot( !l2dist(new_ugdata$udata[new_ugdata$uidx,], new_genodata) ) # 调试用，待会删
      
      ## 连接
      ## 默认使用CS连接法
      if (ncol(new_genodata) <= 30 && !i_missing)
      {
        new_out <- ligation2LCSEM(genodata = new_genodata, udataAB = new_ugdata,
                                  ResA = P_Res[[ 2 * i - 1]],
                                  ResB = P_Res[[ 2 * i ]],
                                  rho = rho, truncthres2 = truncthres2)
      } else {
        new_out <- fastligationEM(ResA = P_Res[[ 2 * i - 1]], ResB = P_Res[[ 2 * i ]], udataAB = new_ugdata, 
                                  truncthres2 = truncthres2)
      }
      New_P_Res[[i]] <- new_out
      New_genodata_list[[i]] <- new_genodata
      New_ugdata_list[[i]] <- new_ugdata
    }
    if (odds) {
      New_P_Res[[new_p_len + 1]] <- last(P_Res)
      New_genodata_list[[new_p_len + 1]]<- last(genodata_list)
      New_ugdata_list[[new_p_len + 1]] <- last(ugdata_list)
      newpconfig[newblocks, ] <- pconfig[blocks,]
    }
    blocks <- newblocks
    pconfig <- newpconfig
    P_Res <- New_P_Res
    genodata_list <- New_genodata_list
    ugdata_list <- New_ugdata_list
  }
  
  out <- New_P_Res[[1]]
  if (inferphase) {
    out$HtDiplo <- as.HtDiplo.HtPhase(out$HtPhase)
  }  else {
    out$HtPhase <- NULL
  }
  
  return(out)
}


# Extra-Long Phasing ------------------------------------------------------
CSHAP_2SNP <- function(genodata, truncthres2 = 1e-3, max_block_size = 10, inferphase = TRUE, ligation_step = 0){
  ## 目前仅仅在每一段上使用了加权，连接的时候还没有加权
  ## 由于新方法只需要使用Phase信息，这里强制推断个体
  ## youligation_step:
  q <- ncol(genodata)
  N_Inds <- nrow(genodata)
  pconfig <- partitionconfig(q, max_block_size = max_block_size)
  
  blocks <- nrow(pconfig)
  genodata_list <- vector(mode = 'list', length = blocks)
  names(genodata_list) <- paste0('block', 1:blocks)
  P_Res <- ugdata_list <- genodata_list
  
  Allhaploi <-  haplo(pconfig[1,3])
  .RM_Mati <- haplo2RM(Allhaploi)
  for (i in seq_len(blocks)){
    genodata_list[[i]] <- genodata[, pconfig[i,1] : pconfig[i,2], drop = FALSE]
    ugdata_list[[i]] <- uniquegenodata(genodata = genodata_list[[i]], MARGIN = 1, detail = TRUE)
    if (i <= blocks - 2) {
      resi <- CSHAPEM(genodata = genodata_list[[i]], Allhaplo = Allhaploi, truncthres2 = truncthres2,
                         udata = ugdata_list[[i]], inferphase = TRUE, .RM_Mat = .RM_Mati)
    } else {
      resi <- CSHAPEM(genodata = genodata_list[[i]], Allhaplo = NULL, truncthres2 = truncthres2,
                         udata = ugdata_list[[i]], inferphase = TRUE, .RM_Mat = NULL)
    }
    if (!ligation_step & is.null(resi$HtDiplo)) resi$HtDiplo <- as.HtDiplo.HtPhase(resi$HtPhase)
    P_Res[[i]]  <- resi
    
    ## 这里需要让每一段都只输出最大概率的那个
    ## 并且应该让h1 <= h2
  }
  
  ### 连接步：将相邻的原子块连接为一个更大的块，防止单个块内的杂合位点太少
  if (ligation_step > 0){
    for (il in 1:ligation_step) {
      last_ligation <- il == ligation_step ## 最后一步时，改用2SNP
      newblocks <- ceiling(blocks / 2)
      odds <- blocks %% 2
      new_p_len <- blocks %/% 2
      newpconfig <- matrix(0, ncol = 3, nrow = newblocks)
      colnames(newpconfig) <- colnames(pconfig)
      
      New_P_Res <- vector('list', newblocks)
      names(New_P_Res) <- 1:newblocks
      New_ugdata_list <- New_P_Res 
      for (i in seq_len(new_p_len)){
        ## 连接 newpconfig
        newpconfig[i, 1] <- pconfig[2*i-1, 1]  
        newpconfig[i, 2] <- pconfig[2*i, 2]
        newpconfig[i, 3] <- pconfig[2*i-1, 3]  + pconfig[2*i, 3]
        
        new_ugdata <- ligation.udata(U1 = ugdata_list[[2 * i - 1]],
                                     U2 = ugdata_list[[2 * i]])
        
        ## 重点：连接相邻块
        new_out <- fastligationEM(ResA = P_Res[[ 2 * i - 1]], ResB = P_Res[[ 2 * i ]], udataAB = new_ugdata, truncthres2 = truncthres2,
                                  uniqueDiplo_in = FALSE, uniqueDiplo_out = last_ligation, transMat = FALSE)
        
        if (last_ligation) new_out$Diplo <- as.Diplo.HtDiplo(new_out$HtDiplo)
        New_P_Res[[i]] <- new_out
        New_ugdata_list[[i]] <- new_ugdata
      }
      
      if (odds) {
        New_P_Res[[new_p_len + 1]] <- last(P_Res)
        New_ugdata_list[[new_p_len + 1]] <- last(ugdata_list)
        newpconfig[newblocks, ] <- pconfig[blocks,]
      }
      
      blocks <- newblocks
      pconfig <- newpconfig
      P_Res <- New_P_Res
      ugdata_list <- New_ugdata_list
      if (blocks <= 2) break
    }
    
    P_Res[[length(P_Res)]]$HtDiplo <- as.HtDiplo.HtPhase(P_Res[[length(P_Res)]]$HtPhase)
  }  
  
  #### debug #####
  # DEBUG <- isFALSE(HtD_pdata)
  
  P2_Res <- vector(mode = 'list', length = blocks - 1)   ## 这里记录的是相邻连接的情况
  names(P2_Res) <- str_c('block', 1:(blocks - 1), '_', 2:blocks)
  ugdata2_list <- P2_Res
  for (i in seq_len(blocks - 1)){
    new_ugdata <- ligation.udata(U1 = ugdata_list[[i]],
                                 U2 = ugdata_list[[i + 1]])
    
    P2_Res[[i]] <- fastligationEM(ResA = P_Res[[i]], ResB = P_Res[[i + 1]], udataAB = new_ugdata, 
                                  uniqueDiplo_in = FALSE, uniqueDiplo_out = TRUE, truncthres2 = truncthres2, transMat = TRUE)
    ugdata2_list[[i]] <- new_ugdata
    P2_Res[[i]]$Diplo <- as.Diplo.HtDiplo(P2_Res[[i]]$HtDiplo)
  }
  
  
  lig_ugdata <- function(iA, iB) {
    ## 将任意两块连接起来得到的udata，使用动态存储结果
    if (iB == iA + 1){
      res <- ugdata2_list[[iA]]
    }else{
      res <- ligation.udata(U1 = ugdata_list[[iA]], U2 = ugdata_list[[iB]])
    }
    return(res)
  }
  
  
  lig_block <- memoise(function(iA, iB){
    ## 将任意两块进行连接的结果，注意这里iA iB不一定相邻
    if (iB == iA + 1){
      res <- P2_Res[[iA]]
    } else {
      res <- fastligationEM(ResA = P_Res[[iA]], ResB = P_Res[[iB]], udataAB = lig_ugdata(iA, iB), 
                           uniqueDiplo_in = FALSE, uniqueDiplo_out = TRUE, 
                           truncthres2 =  truncthres2, transMat = TRUE)
    }
    return(res)
  })
  
  Prob_adj_lig <- function(iA, i, j){
    ## 只计算相邻两块的连接，这里iB = iA + 1
    return(P2_Res[[iA]]$transMatjP[i, j])
  }
  ## 第 A 块 的单体型i 与 第B块的单体型j连接在一起组成新的单体型的概率
  Prob_lig <- memoise(function(iA, iB, i, j) {
    if (iB == iA + 1){
      return(P2_Res[[iA]]$transMatjP[i, j])
    } else {
      ligres <- lig_block(iA, iB)
      return(ligres$transMatjP[i, j])
    }
  })
  
  Odds_adj_lig <- function(iA, i_Ind){
    ## 对于第i个个体，计算其 iA 连接到 iB 块上时，顺式 / 反式 的优势比
    ## 这里顺式定义为 较小的单体型A与较小的单体型B连接，反式相反
    ## 要求A块和B块的结果HtDiplo的每一行是排过序的
    hA <- P_Res[[iA]]$HtDiplo$Diplo[i_Ind,]
    hB <- P_Res[[iA + 1]]$HtDiplo$Diplo[i_Ind,]
    
    Prob_cis <- Prob_adj_lig(iA, i = hA[1], j = hB[1]) * Prob_adj_lig(iA, i = hA[2], j = hB[2])
    Prob_trans <- Prob_adj_lig(iA, i = hA[1], j = hB[2]) * Prob_adj_lig(iA, i = hA[2], j = hB[1])
    ### 注：这里可以做文章，比如使用2SNP的那种算法
    if (Prob_trans == 0) {
      odd_res <- Inf
    } else {
      odd_res <- log(Prob_cis / Prob_trans)
    }
    return(odd_res)
  }
  
  Odds_lig <- memoise(function(iA, iB, i_Ind){
    ## 更加一般的连接，不限于相邻块
    if (iB == iA + 1){
      return(Odds_adj_lig(iA, i_Ind))
    }
    hA <- P_Res[[iA]]$HtDiplo$Diplo[i_Ind,]
    hB <- P_Res[[iB]]$HtDiplo$Diplo[i_Ind,]
    Prob_cis <- Prob_lig(iA, iB, i = hA[1], j = hB[1]) * Prob_lig(iA, iB, i = hA[2], j = hB[2])
    Prob_trans <- Prob_lig(iA, iB, i = hA[1], j = hB[2]) * Prob_lig(iA, iB, i = hA[2], j = hB[1])
    ### 注：这里可以做文章，比如使用2SNP的那种算法
    if (Prob_trans == 0) {
      odd_res <- Inf
    } else {
      odd_res <- log(Prob_cis / Prob_trans)
    }
    # if (abs(odd_res) <= 5){
    #   ## 信息量不够，还需要探索更多位点
    #   genodata[i_Ind, pconfig[iA, 1] : pconfig[iA, 2] ]
    #   genodata[i_Ind, pconfig[iB, 1] : pconfig[iB, 2] ]
    # }
    
    return(odd_res)
  })
  
  ### 由于这里纯合子无信息，会影响后续连接，所以需要对一个个体分析哪些块上是纯合子。
  ##一个 T x m 矩阵，m表示块数，每个元素表示此个体在这一块上是否为1。
  ## 为了防止混淆，以下用i表示个体，j表示块
  Homo_mat <- Matrix(FALSE, ncol = blocks, nrow = N_Inds, dimnames = list(rownames(genodata), names(P_Res)))
  for (j in 1:blocks){
    homo_uidx <- !rowSums(ugdata_list[[j]]$udata == 1)
    homo_uidx[ugdata_list[[j]]$uidx]  ->  Homo_mat[, j]
  }
  Odds_mat <- Matrix(NA, ncol = blocks - 1, nrow = N_Inds, dimnames = list(rownames(genodata), names(P2_Res)))
  
  #### 开始对个体进行推断
  ## 先设计一个矩阵存储最终结果
  HtDiplo_all <- matrix(0, ncol = blocks, nrow = 2 * N_Inds)
  Diplo_all <- Diplo(matrix(0, ncol = q, nrow = 2 * N_Inds))
  rownames(HtDiplo_all) <- str_c(rep(1:N_Inds, each = 2), '_', rep(1:2, N_Inds))
  rownames(Diplo_all) <- rownames(HtDiplo_all)
  HtDiploi <- matrix(0, ncol = blocks, nrow = 2)
  for (i in 1:N_Inds){
    HtDiploi <- matrix(0, ncol = blocks, nrow = 2)
    block_seq <- 1:blocks
    homo_blocks <- block_seq[Homo_mat[i,]]
    hetero_blocks <- block_seq[!Homo_mat[i,]]
    
    if (length(homo_blocks)) {
      HtDiploi[, homo_blocks] <- sapply(homo_blocks, function(o) {
        P_Res[[o]]$HtDiplo$Diplo[i, ]  })
      
      if (length(homo_blocks) == blocks){
        HtDiplo_all[ (2 * i - 1) : (2 * i), ] <- HtDiploi
        next
      }
      
    }
    
    if (length(hetero_blocks) == 1) {
      HtDiploi[, hetero_blocks] <- P_Res[[hetero_blocks]]$HtDiplo$Diplo[i, ]
    } else {
      HtDiploi[, hetero_blocks[1]] <- P_Res[[hetero_blocks[1]]]$HtDiplo$Diplo[i, ]
      cis_j <- 1 
      for (j in seq_along(hetero_blocks[-1])){
        hj <- hetero_blocks[j]
        hj_1 <- hetero_blocks[j + 1]
        oddsij <- Odds_lig(iA = hj, iB = hj_1, i_Ind = i)
        
        Odds_mat[i, hj] <- oddsij
        ## 后续文章：这里有时候Odds不是很大，那么就需要进一步合并
        cis_j <- cis_j * (2 * (oddsij >= 0) - 1)
        if (cis_j >= 0){
          HtDiploi[, hj_1] <- P_Res[[hj_1]]$HtDiplo$Diplo[i, ] # [i, 1:2]
        } else {
          HtDiploi[, hj_1] <- P_Res[[hj_1]]$HtDiplo$Diplo[i, 2:1]
        }
      }
    }
    
    HtDiplo_all[ (2 * i - 1) : (2 * i), ] <- HtDiploi
  }
  
  ## 现在，试图把HtDiplo_all 转换为 Diplo_all
  ## 后续可以改成函数考虑用C写，速度更快
  for (j in 1:blocks){
    Diplo_all[, pconfig[j, 1] : pconfig[j, 2]] <- t(
      P_Res[[j]]$HtDiplo$ht$ht[, HtDiplo_all[,j], drop = FALSE]
    )
  }
  
  # result <- list(Diplo = Diplo_all, runtime = timeuse, P_Res = P_Res, pconfig = pconfig, Odds_mat = Odds_mat, P2_Res = P2_Res)
  result <- list(Diplo = Diplo_all)
  if (q <= 60) result$HtDiplo <- as.HtDiplo.Diplo(X_D = Diplo_all, sort = TRUE)
  
  return(result)
}
