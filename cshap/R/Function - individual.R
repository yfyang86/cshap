
# Base functions ----------------------------------------------------------

tr <- function(X){
  return(sum(diag(X)))
}

Sigmab <- function(miuA, sigmaA, Nind, Pools) {
  q <- length(miuA)
  # eta(i,k) 和 eta(j,l)的协方差
  cov_eta_ikjl <- function(i, k, j, l) {
    tmp0 <- sigmaA[i, j] * sigmaA[k, l] + sigmaA[k, j] * sigmaA[i, l]
    tmp1 <- tmp0 / (4 * Nind ^ 2 * (Pools - 1))
    tmp2 <-
      ((
        miuA[i] * miuA[j] * sigmaA[k, l] + miuA[i] * miuA[l] * sigmaA[k, j] + miuA[k] *
          miuA[l] * sigmaA[i, j] + miuA[k] * miuA[j] * sigmaA[i, l]
      ) / Pools  + tmp0 / Pools ^ 2
      ) / (16 * Nind ^ 4)
    return(tmp1 + tmp2)
  }
  # eta(i,k) 和 omega(l)的协方差
  cov_eta_omega <- function(i, k, l) {
    ans <- (miuA[i] * sigmaA[k, l] + miuA[k] * sigmaA[i, l]) / (8 * Nind ^ 3 *
                                                                  Pools)
    return(ans)
  }
  S22 <- sigmaA / (4 * Nind ^ 2 * Pools)
  S32 <- matrix(0, ncol = q, nrow = q * (q - 1) / 2)
  currentrow32 <- 1
  for (i in 1:(q - 1)) {
    covb3ib2 <- matrix(0, ncol = q, nrow = (q - i))
    for (k in 1:(q - i)) {
      for (l in 1:q) {
        covb3ib2[k, l] <- cov_eta_omega(i, i + k, l)
      }
    }
    S32[currentrow32:(currentrow32 + q - i - 1) ,] <- covb3ib2
    currentrow32 <- currentrow32 + q - i
  }
  S23 <- t(S32)
  S33 <- matrix(0, ncol = q * (q - 1) / 2, nrow = q * (q - 1) / 2)
  currentrow33 <- 1
  for (i in 1:(q - 1)) {
    currentcol33 <- 1
    for (j in 1:(q - 1)) {
      covb3ib3j <- matrix (0, nrow = q - i , ncol = q - j)
      for (k in 1:(q - i)) {
        for (l in 1:(q - j)) {
          covb3ib3j[k, l] <- cov_eta_ikjl(i, i + k, j, j + l)
        }
      }
      S33[currentrow33:(currentrow33 + q - i - 1), currentcol33:(currentcol33 +
                                                                   q - j - 1)] <- covb3ib3j
      currentcol33 <- currentcol33 + q - j
    }
    currentrow33 <- currentrow33 + q - i
  }
  return(rbind(cbind(S22, S23), cbind(S32, S33)))
}

###从Sigma0和omega直接构造输出b_hat误差的模长平方的期望和方差，标准差###
Sigma0b <- function(miu0, sigma0, Nind, Pools) {
  ob <-  Sigmab(miuA = 2 * Nind * miu0,
                sigmaA = 2 * Nind * sigma0,
                Nind,
                Pools)
  varob <- 2 * tr(ob %*% ob)
  return(c(
    ErrorMean = tr(ob),
    ErrorVar = varob,
    ErrorSd = sqrt(varob)
  ))
  
}
## 在一个矩阵中寻找匹配的向量, MARGIN默认为1 （按照行匹配）
match_v_old <- function(x, tab, MARGIN = 1){
  stopifnot(dim(tab)[-MARGIN] == length(x)) ## 长度不匹配
  idx <- apply(tab, MARGIN = MARGIN, function(o){all(o == x)
  })
  return(idx)
}

## 返回匹配的坐标
which_v <- function(x, tab, MARGIN = 1){
  idx <- match_v(x = x, tab = tab, MARGIN = MARGIN)
  return(which(idx))
}

## 返回匹配的元素个数
count_v <- function(x, tab, MARGIN = 1){
  idx <- match_v(x = x, tab = tab, MARGIN = MARGIN)
  return(sum(idx))
}

## 将一组数据转化为unique数据

uniquedata_old <- function(tab, MARGIN = 1){
  udata <- unique(tab, MARGIN = MARGIN)
  count <- apply(udata, MARGIN = MARGIN, count_v, tab = tab)  # 瓶颈
  list(udata = udata, count = count)
}

## 使用RCpp写的版本，速度大幅提升200倍

match_v <- function(x, tab, MARGIN = 1){
  stopifnot(dim(tab)[-MARGIN] == length(x)) ## 长度不匹配
  if (MARGIN == 1){
    idx <- match_row(x, tab)
  } else if (MARGIN == 2){
    idx <- match_col(x, tab)
  }
  
  return(idx)
}




uniquedata <- function(tab, MARGIN = 1, detail = FALSE){
  udata <- unique(tab, MARGIN = MARGIN)
  if (MARGIN == 1)  {
    if (!detail) {
      count <- count_row_mat(udata, tab)
      return(list(udata = udata, count = count))
    } else {
      countres <- count_row_mat_detail(udata, tab)
      N <- nrow(tab)
      uidx <- integer(N)
      for (i in seq_along(countres$count)){
        uidx[countres$matchidx[[i]]] <- i
      }
      
      return(list(
        udata = udata,
        count = countres$count,
        matchidx = countres$matchidx,
        uidx = uidx
      ))
    }
  }  else
    if (MARGIN == 2)  {
      if (!detail) {
        count <- count_col_mat(udata, tab)
        return(list(udata = udata, count = count))
      } else {
        countres <- count_col_mat_detail(udata, tab)
        N <- ncol(tab)
        uidx <- integer(N)
        for (i in seq_along(countres$count)){
          uidx[countres$matchidx[[i]]] <- i
        }
        return(list(
          udata = udata,
          count = countres$count,
          matchidx = countres$matchidx,
          uidx = uidx
        ))
      }
    }
  
}

### 仅当 q < 19 生效，要求行名或者列名是唯一的。

## 注意：要求genodata 对应margin的名字是唯一的
uniquegenodata <- function(genodata, MARGIN = 1, detail = FALSE, missing = FALSE){
  if (MARGIN == 1) {
    q <- ncol(genodata)
    if (q > 19) return (uniquedata(tab = genodata, MARGIN = MARGIN, detail = detail))
    Nind <- nrow(genodata)
    if (is.null(rownames(genodata)))
      rownames(genodata) <- seq_len(Nind)
  } else {
    q <- nrow(genodata)
    if (q > 19) return (uniquedata(tab = genodata, MARGIN = MARGIN, detail = detail))
    Nind <- ncol(genodata)
    if (is.null(colnames(genodata)))
      colnames(genodata) <- seq_len(Nind)
  }
  
  base_g <- (3 + missing)^((q-1) : 0)
  ## 由于缺失位点是用3表示的。
  
  if (MARGIN == 1) {
    ndata <- genodata %*% base_g
    rownames(ndata) <- 1:nrow(genodata)
  } else {
    ndata <- as.matrix(colSums(genodata * base_g))
    rownames(ndata) <- 1:ncol(genodata)
  }
  
  undata <- unique(ndata)
  
  if (MARGIN == 1) {
    udata <- genodata[as.integer(rownames(undata)) ,, drop = FALSE]
  } else {
    udata <- genodata[, as.integer(rownames(undata)), drop = FALSE]
  }
  
  count <- countint(undata, ndata)
  outres <- list(udata = udata, count = count)
  # if (!detail){
  # return(outres)
  # } else {
  #   # countres <- countint_detail(undata, ndata)
  #   uidx <- match(ndata, undata)
  #   outres$uidx <- uidx
  #   # return(list(udata = udata, count = countres$count, matchidx = countres$matchidx))
  # }
  if (detail) {
    uidx <- match(ndata, undata)
    outres$uidx <- uidx
  }
  return(outres)
}



# 限制 q < 30
uniquehaplodata <- function(genodata, MARGIN = 1, detail = FALSE){
  if (MARGIN == 1) {
    q <- ncol(genodata)
    if (q > 29) return (uniquedata(tab = genodata, MARGIN = MARGIN, detail = detail))
    Nind <- nrow(genodata)
    if (is.null(rownames(genodata)))
      rownames(genodata) <- seq_len(Nind)
  } else {
    q <- nrow(genodata)
    if (q > 29) return (uniquedata(tab = genodata, MARGIN = MARGIN, detail = detail))
    Nind <- ncol(genodata)
    if (is.null(colnames(genodata)))
      colnames(genodata) <- seq_len(Nind)
  }
  
  base2 <- 2^( (q-1) : 0)
  
  if (MARGIN == 1) {
    ndata <- genodata %*% base2
  } else {
    ndata <- as.matrix(colSums(genodata * base2))
  }
  
  undata <- unique(ndata)
  
  if (MARGIN == 1) {
    udata <- genodata[rownames(undata) ,]
  } else {
    udata <- genodata[, rownames(undata)]
  }
  
  outres <- list(udata = udata, count = count)
  if (detail){
    uidx <- match(ndata, undata)
    outres$uidx <- uidx
  }
  return(outres)
  
}


boot_rho <- function(genodata, Nind = 1, nboot = 100, check = TRUE){
  tmp <- bootstrap(1:nrow(genodata), nboot , function(o) {
    estimate.rho(genodata[o, ], Nind = Nind, check = check)
  })
  rho <- estimate.rho(genodata, Nind = Nind, check = check)
  rho_sd <- sd(tmp$thetastar)
  rho_L <- (abs(rho)  - 1.645 * rho_sd) > 0 
  # rho_Z <- rho / rho_sd
  return(rho_L)
}

estimate.rho <- function(genodata, Nind, check = TRUE){
  if (check) genodata[genodata == 3] <- NA
  VarG <- apply(genodata,2,var, na.rm = TRUE)
  MeanG <- colMeans(genodata, na.rm = TRUE)
  e.rho <- (sum(VarG)/sum(MeanG*(1 - MeanG/(2*Nind)))) -1
  return(e.rho)
}

estimate_rho <- function(genodata, N = 1){
  return(estimate.rho(genodata = genodata, Nind = N, check = TRUE))
}

estimate.MAF <- function(genodata, Nind, rho = 0, check = TRUE) {
  ## check = TRUE 意味着需要把3换成NA
  if (check) genodata[genodata == 3] <- NA
  SigmaG <- var(genodata, na.rm = TRUE)
  VarG <- diag(SigmaG)
  MeanG <- as.matrix(colMeans(genodata, na.rm = TRUE))
  omega_hat <- MeanG / (2 * Nind)
  Sigma_hat <- SigmaG / (2 * Nind)
  colnames(MeanG) <- 'freq'
  eta_rho <- Sigma_hat / (1 + rho) + omega_hat %*% t(omega_hat)
  return(list(
    omega_hat = omega_hat,
    eta_hat = eta_rho,
    HWE = abs(rho) <= 0.001,
    rho = rho
  ))
}


estimate.b <- function(genodata, e.MAF = NULL,...){
  if (is.null(e.MAF)) {e.MAF <- estimate.MAF(genodata = genodata,...) }
  omega <- e.MAF$omega_hat
  eta <- e.MAF$eta_hat
  q <- length(omega)
  b1 <- 1
  b2 <- omega
  b3 <- matrix(0,ncol = 1,nrow = q*(q-1)/2)
  b3currentrow <- 0 
  rownames(b3) <- 1:(q*(q-1)/2)
  colnames(b3) <- 'freq'
  for (i in 1:(q-1)){
    for (j in (i+1):q){
      b3currentrow <- b3currentrow + 1
      b3[b3currentrow,] <- eta[i,j]
      rownames(b3)[b3currentrow] <- paste0('loci',i,'_',j)
    }
  }
  b <- rbind(b1,b2,b3)
  return(b)
}



####衡量一次实验结果的指标：haplotype distance (p0,estimate.p) 单体型分布估计和真实值之间的平均欧氏距离 ###
htdist <- function(p0,estimate.p){
  if (class(p0) != 'Sparsep')   p0 <- as.Sparsep(p0)
  if (class(estimate.p)!= 'Sparsep')  estimate.p <- as.Sparsep(estimate.p)
  q <- p0$q
  r <- 2^q
  pdist <- sqrt(sum((p0$p - estimate.p$p)^2)/r)
  su1 <- unique(p0$su_p)
  su2 <- unique(estimate.p$su_p)
  TP <- sum(su1 %in% su2)
  FN <- length(su1) - TP
  FP <- length(su2) - TP
  TN <- r -length(su1) -length(su2) + TP
  Precision <- TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F1 <- 2*TP / (2*TP + FP + FN)
  #### Confusion Matrix
  Conf.Mat <-  cbind(c(TP,FP),c(FN,TN))
  colnames(Conf.Mat) <- c('Actual_1','Actual_0')
  rownames(Conf.Mat) <- c('Predict_1','Predict_0')
  ### 有效累积概率 EAP
  EAP <- sum(estimate.p$p[intersect(su1,su2)])
  result <- list(q=q,su_p0=su1,su_estimate.p=su2,Conf.Mat=Conf.Mat,pdist=pdist,Precision=Precision,Recall=Recall,F1=F1,EAP=EAP)  
  return(result)
}




# EM functions ------------------------------------------------------------
E_Step0 <- function(genotype){
  q <- length(genotype)
  ho2 <- genotype == 2
  ho0 <- genotype == 0
  ho1 <- genotype == 1
  si <- sum(ho1)
  if (si <= 1) {
    result <- as.integer(as.htn.SNPvt(rbind(floor(genotype / 2), ceiling(genotype / 2)),
                                      q = q))
    result <- matrix(c(result, Prob = 1) , nrow = 1)
    return(result)
  }
  
  base <- 2^((q-1):0)
  base1 <- base[ho1]
  excluded <- sum(base[ho2])
  base1_2 <- base1[-1]
  if (si == 2){
    htnsettmp <- colSums(haplo(si) * base1) + excluded + 1
    h1 <- htnsettmp[1:2]
    h2 <- htnsettmp[4:3]
  } else {
    h1 <- colSums(haplo(si - 1) * base1_2) + excluded + 1
    h2 <- rev(h1) + base1[1]
  }
  Prob <- 1/ 2^(si-1)
  # colSums(haplo(si) * base1) + excluded + 1
  # genotype1 <- genotype[ho1] of course, 1 1 1
  diploresult <- cbind(h1, h2, Prob)
  return(diploresult)
}

E_Step <- function(genotype, freq, Allhaplo = NULL) {
  ## 这里的E_Step 其实是第一步E_Step ，分工明确。不进行截断
  q <- length(genotype)
  ho2 <- genotype == 2
  ho0 <- genotype == 0
  ho1 <- genotype == 1
  si <- sum(ho1)
  if (si <= 1) {
    result <- as.integer(as.htn.SNPvt(rbind(floor(genotype / 2), ceiling(genotype / 2))))
    result <- matrix(c(result, 1) , nrow = 1)
    return(result)
  }
  base <- 2^((q-1):0)
  excluded <- sum(base[ho2])
  if (is.null(Allhaplo)) Allhaplo <- haplo(q)
  # htnset <- sort(supp(freq)) 无需 sort
  htnset <- which(freq > 0)
  haplo_set <- Allhaplo[, htnset, drop = FALSE]
  ## h1idx <- apply(haplo_set, 2, function(o){all(o[ho2] == 1) & all(o[ho0] ==0)}) 
  h1idx <- (colSums(haplo_set[ho2,, drop = FALSE]) == sum(ho2)) & (colSums(haplo_set[ho0,, drop = FALSE]) == 0) ## Faster
  h1 <- htnset[h1idx]
  if (!any(h1idx)) return(E_Step0(genotype = genotype))
  base1 <- base[ho1]
  
  freq1 <- freq[h1]
  # genotype1 <- genotype[ho1]
  # genotype1 <- rep(1, si)
  haploset1 <- Allhaplo[ho1, h1, drop = FALSE]
  Prob <- h2 <- integer(length(h1))
  # diplocompati <- h2 ## 是否相容
  # compati 记录了相容情况——2:相容，1:次相容，0:不相容
  # h2 ## 概率
  # compati <- h2
  for (i in seq_along(h1)) {
    htni <- h1[i]
    if (any(htni == h2)) next # Unique Config !!!
    reshaplo <- 1 - haploset1[,i]
    if (any(reshaplo < 0 | reshaplo >= 2)) {
      # compati[i] <- 0
      next
    } 
    reshtn <- sum(base1 * reshaplo) + excluded + 1
    h2[i] <- reshtn
    Prob[i] <- freq1[i] * (2 * freq[reshtn] + min(freq1[i]) * (!freq[reshtn]) / 5 ) ##  经过测试，取10比较好
    # diplocompati[i] <- reshtn %in% htnset
  }
  sumprob <- sum(Prob)
  if (!sumprob) return(E_Step0(genotype = genotype)) ## 鲁棒
  Prob <- Prob / sumprob
  diploresult <- cbind(h1, h2, Prob)[Prob >0 ,,drop = FALSE]
  return(diploresult)
}

trunc_E_Result <- function(E_Result, Prob = E_Result[,3], truncthres2 = 1e-3){
  # Prob <- E_Result[, 3]
  Prob_Tidx <- Prob < truncthres2
  if (all(Prob_Tidx)) {
    return(E_Result)
  } else {
    Prob[Prob_Tidx] <- 0
    Prob <- Prob / sum(Prob)
    E_Result[, 3] <- Prob
    E_Result <- E_Result[!Prob_Tidx , , drop = FALSE]
    
    if (nrow(E_Result) == 1){
      E_Result[,1:2] <- sort(E_Result[,1:2])
    }
    
    return(E_Result)
  }
  
}


E_Step_Miss0 <- function(genotype){
  q <- length(genotype)
  ho2 <- genotype == 2
  ho0 <- genotype == 0
  ho1 <- genotype == 1
  ho3 <- genotype == 3
  si <- sum(ho1)
  base <- 2^((q-1):0)
  if (si <= 1) {
    genotype2 <- genotype
    genotype2[ho3] <- 0
    result <- as.integer(as.htn.SNPvt(rbind(floor(genotype2 / 2), ceiling(genotype2 / 2)),
                                      q = q))
    htnsresult <- matrix(result , nrow = 1)
  } else {
    base1 <- base[ho1]
    excluded <- sum(base[ho2])
    base1_2 <- base1[-1]
    if (si == 2){
      htnsettmp <- colSums(haplo(si) * base1) + excluded + 1
      h1 <- htnsettmp[1:2]
      h2 <- htnsettmp[4:3]
    } else {
      h1 <- colSums(haplo(si - 1) * base1_2) + excluded + 1
      h2 <- rev(h1) + base1[1]
    }
    htnsresult <- cbind(h1, h2)
  }
  
  ### 处理缺失位点
  base3 <- base[ho3]
  for (i in length(base3):1){
    htnsresult <- add_miss_loci(htnsresult, miss_base = base3[i])
  }
  
  Prob <- 1/ (1 + (htnsresult[,1] == htnsresult[,2] ))
  diploresult <- cbind(htnsresult, Prob/sum(Prob))
  
  return(diploresult)
}


E_Steps_mixed0 <- function(genodata, miss_row){
  ### 混合了缺失与没缺失的E步
  E_Results <- vector('list', length = nrow(genodata))
  if(any(!miss_row)) E_Results[!miss_row] <- alply(.data = genodata[!miss_row, , drop = FALSE], .margins = 1, .fun = E_Step0, .dims = FALSE)
  if(any(miss_row)) E_Results[miss_row] <- alply(.data = genodata[miss_row, , drop = FALSE], .margins = 1, .fun = E_Step_Miss0, .dims = FALSE)
  return(E_Results)
}

## 从一个双体型配置中加入一个缺失位点
sub_add_miss_loci <- function(htns, miss_base){
  ## 由于区分纯合子与杂合子比较麻烦，干脆不予区分，纯合子的第三行用0来填充
  if (htns[1] == htns[2]){
    ## 纯合子只有2种
    res <- matrix(htns, byrow = TRUE, nrow = 4, ncol = 2, dimnames = list(c(), c('h1', 'h2')))
    res[2,] <- res[2,] + c(0, miss_base)
    res[3,] <- 0
    res[4,] <- res[4,] + miss_base
  } else {
    res <- matrix(htns, byrow = TRUE, nrow = 4, ncol = 2, dimnames = list(c(), c('h1', 'h2')))
    res[2,] <- res[2,] + c(0, miss_base)
    res[3,] <- res[3,] + c(miss_base, 0)
    res[4,] <- res[4,] + miss_base
  } 
  return(res)
}

## 从一组双体型配置中加入一个缺失位点
add_miss_loci <- function(htnsresult, miss_base){
  n <- nrow(htnsresult)
  Res <- matrix(0, nrow = 4*n ,ncol = 2,  dimnames = list(c(), c('h1', 'h2')))  
  for (i in 1:n){
    sub_add_miss_loci(htnsresult[i,], miss_base) -> Res[(4*i - 3) : (4*i), ]
  }
  return(Res[Res[,1] > 0, , drop = FALSE])
}


### 考虑缺失位点的M_Step打算用C写，速度应该有明显提升
M_Step_miss <- function(E_Result, r){
  newfreq <- numeric(r)
  if (E_Result[1,1] == E_Result[1,2] && nrow(E_Result) == 1) {
    newfreq[E_Result[1,1]] <- 1
    return(newfreq)
  } 
  
  for (i in 1:nrow(E_Result)){
    
    prob <- E_Result[i, 3] / 2## Prob 应该已经被归一化了
    newfreq[E_Result[i, 1]] <- newfreq[E_Result[i, 1]] + prob
    newfreq[E_Result[i, 2]] <- newfreq[E_Result[i, 2]] + prob
    
  }
  # newfreq[c(E_Result[,1:2])] <- E_Result[,3] / 2
  return(newfreq)
}

M_Steps_mixed <- function(E_Results, r, miss_row){
  Res <- matrix(0, ncol = length(E_Results), nrow = r)
  if(any(!miss_row)) Res[,!miss_row] <- sapply(E_Results[!miss_row], M_Step, r = r)
  if(any(miss_row)) Res[,miss_row] <- sapply(E_Results[miss_row], M_Step_miss, r = r)
  return(Res)
}

M_Steps <- function(E_Results, r){
  res <- sapply(E_Results, M_Step, r = r)
  return(res)
}

E_Step_Miss <- function(genotype, freq, Allhaplo = NULL) {
  ## 这里的E_Step 其实是第一步E_Step ，分工明确。不进行截断
  q <- length(genotype)
  ho2 <- genotype == 2
  ho0 <- genotype == 0
  ho1 <- genotype == 1
  ho3 <- genotype == 3
  # si <- sum(ho1)
  base <- 2^((q-1):0)
  base1 <- base[ho1]
  excluded <- sum(base[ho2])
  if (is.null(Allhaplo)) Allhaplo <- haplo(q)
  # htnset <- sort(supp(freq)) 无需 sort
  htnset <- which(freq > 0)
  haplo_set <- Allhaplo[, htnset, drop = FALSE]
  ## h1idx <- apply(haplo_set, 2, function(o){all(o[ho2] == 1) & all(o[ho0] ==0)}) 
  h1idx <- (colSums(haplo_set[ho2,, drop = FALSE]) == sum(ho2)) & (colSums(haplo_set[ho0,, drop = FALSE]) == 0) ## Faster
  h1 <- htnset[h1idx]
  if (!any(h1idx)) return(E_Step_Miss0(genotype = genotype))
  if (!any(ho1)){
    ## h1 中的元素两两组合，但是纯合子只算一种
    n_h1 <- sum(h1idx)
    htnsresult <- matrix(0, nrow = n_h1 * (n_h1 + 1) / 2, ncol = 2, dimnames = list(c(), c('h1', 'h2')))
    crow <- 1
    for (i in 1:n_h1){
      for (j in i:n_h1){
        htnsresult[crow, ] <- h1[c(i, j)]
        crow <- crow + 1
      }      
    }
    htnsresult[,1] -> h1
    htnsresult[,2] -> h2
    Prob <- freq[h1] * (freq[h2] + min(freq[h1]) * !freq[h2] / 5 )  / (1 + (h1 == h2))
    
  } else {
    haploset1 <- Allhaplo[ho1, h1, drop = FALSE]
    h2 <- integer(length(h1))
    for (i in seq_along(h1)) {
      htni <- h1[i]
      if (any(htni == h2)) next # Unique Config !!!
      reshaplo <- 1 - haploset1[,i]
      if (any(reshaplo < 0 | reshaplo >= 2)) {
        # compati[i] <- 0
        next
      } 
      reshtn <- sum(base1 * reshaplo) + excluded + 1
      h2[i] <- reshtn
    }
    
    h1 <- h1[h2 > 0]
    h2 <- h2[h2 > 0]
    htnsresult0 <- cbind(h1, h2) ## 还没有考虑缺失位点时候暂时的htn结果
    
    ### 由于h1给定了在supp中取，只需要考虑h2缺失位点的补全即可
    
    base3 <- base[ho3]
    miss_base <-  base3 %*% haplo(sum(ho3))
    n_htn_0 <- nrow(htnsresult0)
    htnsresult <- matrix(0, nrow = length(miss_base) * n_htn_0, ncol = 2, dimnames = list(c(), c('h1', 'h2')))
    htnsresult[, 1] <- h1
    htnsresult[, 2] <- h2
    for (i in seq_along(miss_base)){
      htnsresult[((i - 1)* n_htn_0 + 1): (i * n_htn_0), 2] <- htnsresult[((i - 1)* n_htn_0 + 1): (i * n_htn_0), 2] + miss_base[i]
    }
    htnsresult[,1] -> h1
    htnsresult[,2] -> h2
    Prob <- freq[h1] * (freq[h2] + min(freq[h1]) * !freq[h2] / 5 )  
  }
  
  sumprob <- sum(Prob)
  if (!sumprob) return(E_Step_Miss0(genotype = genotype)) ## 鲁棒
  Prob <- Prob / sumprob
  diploresult <- cbind(h1, h2, Prob)[Prob >0 ,,drop = FALSE]
  return(diploresult)
}


E_Steps_mixed <- function(genodata, freq, miss_row, Allhaplo = NULL){
  ### 混合了缺失与没缺失的E步
  E_Results <- vector('list', length = nrow(genodata))
  
  if(any(!miss_row)) E_Results[!miss_row] <-alply(.data = genodata[!miss_row, , drop = FALSE], .margins = 1,
                                                  .fun = E_Step, freq = freq, Allhaplo = Allhaplo, .dims = FALSE)
  if(any(miss_row)) E_Results[miss_row] <- alply(.data = genodata[miss_row, , drop = FALSE], .margins = 1,
                                                 .fun = E_Step_Miss,
                                                 freq = freq, Allhaplo = Allhaplo, .dims = FALSE)
  return(E_Results)
}

sort_E_Result <- function(E_Result){
  if (nrow(E_Result) ==  1) return(E_Result)
  E_Result[, 3] -> Prob
  E_Result2 <-  aaply(E_Result[,1:2, drop = FALSE], .margins = 1, function(o){sort(o)})
  E_Result[, 1:2] <- E_Result2
  
  ures <- uniquedata(E_Result2, MARGIN = 1, detail = TRUE)
  if (any(ures$count > 1)){
    uProb <- tapply(Prob, ures$uidx, sum)
    E_Result <- cbind(ures$udata, uProb)
  }
  rownames(E_Result) <- NULL
  return(E_Result)
}

E_Step_withdiplo <- function(E_Result, freq, truncthres2 = 1e-3)  {
  ### 若已知可能的双体型配置，则不需要重新推断了
  if(nrow(E_Result) ==  1) 
  {
    E_Result[,3] <- 1
    
    E_Result[,1:2] <- sort(E_Result[,1:2])
    
    return(E_Result)
  }
  Prob <- freq[E_Result[,1]] * freq[E_Result[,2]] / (1 + (E_Result[,1] == E_Result[,2]))
  Prob <- Prob / sum(Prob)
  E_Result[,3] <- Prob
  res <- trunc_E_Result(E_Result = E_Result, Prob = Prob, truncthres2 = truncthres2)
  
  return(res)
}


E_Steps_withdiplo <- function(E_Results, freq, truncthres2 = 1e-3 ){
  res <- llply(.data = E_Results, .fun = E_Step_withdiplo,  freq = freq , truncthres2 = truncthres2)
  return(res)
}

M_Step <- function(E_Result, r){
  newfreq <- numeric(r)
  if (E_Result[1,1] == E_Result[1,2]) {
    newfreq[E_Result[1,1]] <- 1
    return(newfreq)
  } 
  prob <- E_Result[,3] / 2## Prob 应该已经被归一化了
  newfreq[E_Result[,1]] <- prob
  newfreq[E_Result[,2]] <- prob
  # newfreq[c(E_Result[,1:2])] <- E_Result[,3] / 2
  return(newfreq)
}

EM_ind <- function(genotype, p_initial = NULL, Allhaplo = NULL,  truncthres2 = 1e-4, inferphase = FALSE){
  q <- length(genotype)
  r <- 2^q
  genotype[is.na(genotype)] <- 3
  miss <- any(genotype == 3)
  if (is.null(Allhaplo)) Allhaplo <- haplo(q)
  if (miss){
    iE_Step <- E_Step_Miss
    iE_Step0 <- E_Step_Miss0
    iM_Step <- M_Step_miss
  } else {
    iE_Step <- E_Step
    iE_Step0 <- E_Step0
    iM_Step <- M_Step
  }
  if(is.null(p_initial) || is.na(p_initial)){
    freq <- 0
    E_Result <- iE_Step0(genotype = genotype)
  } else {
    freq <- p_initial
    E_Result <- iE_Step(genotype = genotype, freq = freq, Allhaplo = Allhaplo) 
  }
  newfreq <- iM_Step(E_Result, r)
  
  delta <- sum(abs(newfreq - freq))
  # step <- 1
  while (delta > 1e-06){
    freq <- newfreq
    E_Result <- E_Step_withdiplo(E_Result, freq = freq, truncthres2 = truncthres2)
    newfreq <- iM_Step(E_Result, r)
    delta <- sum(abs(newfreq - freq))
  }
  # return(step)
  S_idx <- newfreq > 0
  ht <- Allhaplo[, S_idx, drop = FALSE]
  hf <- haplofreq0(q = q, ht = ht, p = newfreq[S_idx])
  out <- list(freq = hf)
  if (inferphase){
    X_SpPhase <- list(E_Result)
    out$SpPhase <- X_SpPhase
  }
  return(out)
}

fEM <- function(genodata, p_initial = NULL, Allhaplo = NULL, 
                        truncthres2 = 1e-3, 
                        inferphase = FALSE, udata = NULL){
  ## 可以处理缺失位点
  genodata[is.na(genodata)] <- 3
  q <- ncol(genodata)
  r <- 2^q
  N_Inds <- nrow(genodata)
  if (is.null(Allhaplo)) Allhaplo <- haplo(q)
  
  miss_row <- as.logical(rowSums(genodata == 3))
  # if (!any(miss_row)) return(fastEM(genodata = genodata, p_initial = p_initial, Allhaplo = Allhaplo,
  #                                   truncthres2 = truncthres2,  inferphase = inferphase,  udata))
  if (is.null(udata)) udata <- uniquegenodata(genodata, detail = inferphase, missing = any(miss_row))
  
  # 处理缺失的部分 ------
  ugdata <- udata$udata
  ucount <- udata$count
  
  miss_row <- as.logical(rowSums(ugdata == 3))
  
  if (length(ucount) == 1) {
    EM_out <- EM_ind(genotype = ugdata, p_initial = p_initial, Allhaplo = Allhaplo, truncthres2 = truncthres2, inferphase = inferphase)
    if (inferphase) {
      X_HtP <- as.HtPhase.SpPhase(X_SP = EM_out$SpPhase, q = q)
      X_HtP$Phase <- replicate(N_Inds, X_HtP$Phase[[1]], simplify = FALSE)
      EM_out$HtPhase <- X_HtP
    }
    return(EM_out)
  }
  
  if (is.null(p_initial) || is.na(p_initial)) {
    UE_Results <- E_Steps_mixed0(genodata = ugdata, miss_row = miss_row)
    freq <- 0
    # freq <- rep(1 / r, r)
  } else {
    freq <- p_initial
    ## 先把那些可以确定的单体型挑出来——这是非常重要的信息
    confirm.idx  <- (rowSums(ugdata == 1) <= 1)  & !miss_row
    if (any(confirm.idx)){
      confirmeddata <- ugdata[confirm.idx,, drop = FALSE]
      h1s <- as.htn.SNPvt(floor(confirmeddata/2), q = q)
      h2s <- as.htn.SNPvt(ceiling(confirmeddata/2), q = q)
      confirmed.htn <-  unique(c(h1s, h2s))
      freq[confirmed.htn] <- pmax(freq[confirmed.htn], 0.5 / sum(freq > 0), 1/N_Inds) ## 送分的频率要确保存在
      
      E_Res_confirmed <- cbind(h1 = h1s, h2 = h2s, Prob = 1L)
      UE_Results <- vector('list', length = nrow(ugdata))
      UE_Results[confirm.idx] <- lapply(1:sum(confirm.idx), function(o){E_Res_confirmed[o,,drop = FALSE]})
      
      if (any(!confirm.idx))  UE_Results[!confirm.idx] <- E_Steps_mixed(genodata = ugdata[!confirm.idx,, drop = FALSE], freq = freq,
                                                                        miss_row = miss_row[!confirm.idx], Allhaplo = Allhaplo)
    } else {
      UE_Results <- E_Steps_mixed(genodata = ugdata, freq = freq, miss_row = miss_row, Allhaplo = Allhaplo)
    }
  }
  # UE_Results[miss_row]
  
  UM_Res <- M_Steps_mixed(UE_Results, r = r, miss_row = miss_row)
  newfreq <- c(UM_Res %*% ucount) / N_Inds
  delta <- sum(abs(newfreq - freq))
  ### 仅仅考虑那些需要更新的个体
  # genodata_1 <- genodata[indidx, ,drop = FALSE]
  
  interestidx <- which(sapply(UE_Results, nrow) > 1)
  n_interest <- length(interestidx)
  UM_Res_fixed <- UM_Res[, -interestidx, drop = FALSE] %*% ucount[-interestidx]
  # istep <- 1
  UE_Results[miss_row] <- llply(.data = UE_Results[miss_row], .fun = sort_E_Result)
  while (delta > 1e-06 && n_interest){
    # cat(delta, '\n')
    freq <- newfreq
    UE_Results_interest <- UE_Results[interestidx]
    UE_Results_interest <- E_Steps_withdiplo(E_Results = UE_Results_interest, freq = freq, truncthres2 = truncthres2)
    UE_Results[interestidx] <- UE_Results_interest
    new_int_id <- sapply(UE_Results_interest, nrow) > 1 
    new_n_interest <- sum(new_int_id)
    UM_Res_interest <- M_Steps_mixed(UE_Results_interest, r = r, miss_row = miss_row[interestidx])
    if (new_n_interest < n_interest) {
      UM_Res_fixed <- UM_Res_fixed + 
        UM_Res_interest[,!new_int_id, drop = FALSE] %*% ucount[interestidx][!new_int_id]
    }
    newfreq <- c(UM_Res_interest[,new_int_id, drop = FALSE] %*% ucount[interestidx][new_int_id] + 
                   UM_Res_fixed) / N_Inds
    delta <- sum(abs(newfreq - freq))
    interestidx <- interestidx[new_int_id] 
    n_interest <- new_n_interest
    # step <- step + 1
    # deltas <- append(deltas, delta)  仅作为debug功能
  }
  
  newfreq <- newfreq / sum(newfreq)
  S_idx <- newfreq > 0
  ht <- Allhaplo[, S_idx, drop = FALSE]
  hf <- haplofreq0(q = q, ht = ht, p = newfreq[S_idx])
  out <- list(freq = hf)
  
  if (inferphase){
    NewE_Results <- UE_Results[udata$uidx]
    names(NewE_Results) <- rownames(genodata)
    # out$SpPhase <- NewE_Results 
    X_HtPhase <- as.HtPhase.SpPhase2(X_SP = NewE_Results, Sp_Sidx = which(S_idx), ht = ht, q = q)
    # X_HtD <- as.HtDiplo.HtPhase(X_HtPhase)
    # EM_out$HtPhase <- X_HtPhase
    X_HtPhase$ht <- hf
    out$HtPhase <- X_HtPhase
    # out$SpDiplo <- as.SpDiplo.SpPhase(X_SP = NewE_Results, q = q)
  }
  return(out)
  # return(list(freq = newfreq, deltas = deltas))
}

CSHAP <- function(genodata, Nind = 1, 
                  quantile_thres = 0.90, rho = 0, HWE = TRUE, 
                  Allhaplo = NULL, .RM_Mat = NULL){
  Pools <- nrow(genodata)
  q <- ncol(genodata)
  genodata[genodata == 3] <- NA
  if (is.null(Allhaplo)) {Allhaplo = haplo(q)}
  if (is.null(.RM_Mat)) .RM_Mat <- haplo2RM(Allhaplo)
  if (!HWE && !rho){
    ## Adaptive
    if (boot_rho(genodata, Nind = Nind, nboot = 100, check = FALSE)){
      rho <- estimate.rho(genodata, Nind = Nind, check = FALSE)
    } else {
      rho <- 0
    }
  }
  
  tmpe_e_MAF <- estimate.MAF(genodata, Nind = Nind, rho = rho, check = FALSE)
  e_b <- estimate.b(genodata, e.MAF = tmpe_e_MAF, Nind = Nind)
  ### 0 - solution
  restmp0 <- spams.lasso(X = e_b, D = .RM_Mat$phi, return_reg_path = FALSE,
                         lambda1 = 0 , mode = 'L2ERROR', pos = TRUE)
  if (quantile_thres <= 0.009) {
    tmpresult <- restmp0/.RM_Mat$tau
    return(tmpresult/sum(tmpresult))
  }
  tmpMAF <- tmpe_e_MAF
  eta0 <- tmpMAF$eta
  miu0 <- tmpMAF$omega
  sigma0 <- eta0 - miu0%*%t(miu0)
  b0 <- e_b
  ## ju Xingxi  
  ## jiashe Y=\lambda X, X dist Chisquare[k]
  ## es[1]: Mean; es[2]:Var.
  
  es <- Sigma0b(miu0 = miu0,sigma0 = sigma0,Nind = Nind,Pools = Pools)
  estimate_chisquare_df <- 2*es[1]^2/es[2]
  estimate_chisquare_lambda <- 1/2*es[2]/es[1]
  ### Quantile of L2 ERROR SQUARE 
  qerrorsquare <- function(p = 0.95){
    return(estimate_chisquare_lambda * qchisq(p = p,df = estimate_chisquare_df))
  }
  if (Nind == 1) {
    genodata_t <- t(genodata)
    compatibletmp <-
      memoise(
        partial(
          compatible_tdata0,
          genodata = genodata_t,
          Allhaplo = Allhaplo,
          q = q,
          pools = Pools
        )
      )
  } else {
    restmp1 <- spams.lasso(X = e_b,D = .RM_Mat$phi, return_reg_path = FALSE,
                           lambda1 = qerrorsquare(quantile_thres) ,mode='L2ERROR',pos = TRUE)
    tmpresult <- restmp1/.RM_Mat$tau
    return(tmpresult/sum(tmpresult))
  }
  x_min <- sum((.RM_Mat$phi %*% restmp0 - e_b) ^2)
  quantile_seq <- c(seq(0, 0.09, 0.01), seq(0.1, 0.28, 0.02), seq(0.3, 0.96, 0.03))
  quantile_seq <- quantile_seq[quantile_seq <= quantile_thres]
  x_seq <- qerrorsquare(quantile_seq)
  l2errorlist <- filter(tbl_df(data.frame(quantile = quantile_seq, x = x_seq)),x > x_min)
  l2errorlist <- rbind(c(0, x_min), l2errorlist)
  if (nrow(l2errorlist) == 1) {
    tmpresult <- restmp0/.RM_Mat$tau
    return(tmpresult/sum(tmpresult))
  }
  li <- nrow(l2errorlist)
  l2errorlist <- mutate(l2errorlist, li = 1:li)
  # Compatibility Check
  rhapresulttmp <- data.frame(xi = l2errorlist$li, l2errorlist$x, compat = 0)
  rhapresulttmp$l0 <- length(restmp0@i)
  rhaptmpres <- Matrix(0, ncol = 2^q , nrow = nrow(rhapresulttmp))
  if (Nind == 1) { 
    # rhapresulttmp$compat <- sapply(rhapresultallhtn, compose(sum,compatibletmp)) 
    for (i in seq_along(rhapresulttmp$compat)){
      tmpresi <- spams.lasso(X = e_b,D = .RM_Mat$phi,return_reg_path = FALSE,
                             lambda1 = l2errorlist$x[i],mode='L2ERROR',pos = TRUE)
      rhaptmpres[i,] <- tmpresi
      htni <- tmpresi@i + 1
      icomp <- sum(compatibletmp(htni))
      rhapresulttmp$compat[i] <- icomp
      rhapresulttmp$l0[i] <- length(htni)
      if (icomp < rhapresulttmp$compat[max(1,i - 1)]) break
    }
  }
  
  rhapoptresult <- filter(filter(rhapresulttmp, compat == max(compat)), l0 == min(l0))
  tmpresult <- rhaptmpres[max(rhapoptresult$xi),]/.RM_Mat$tau
  return(tmpresult/sum(tmpresult))
}


# CSHAP for individuals ---------------------------------------------------

CSHAP0 <- function(genodata, Nind = 1, rho = 0, HWE = TRUE, Allhaplo = NULL, .RM_Mat = NULL){
  Pools <- nrow(genodata)
  q <- ncol(genodata)
  if (is.null(Allhaplo)) Allhaplo <- haplo(q)
  genodata[genodata == 3] <- NA
    if (!HWE && !rho){
    ## Adaptive
    if (boot_rho(genodata, Nind = Nind, nboot = 100, check = FALSE)){
      rho <- estimate.rho(genodata, Nind = Nind, check = FALSE)
    } else {
      rho <- 0
    }
  }
  tmpe_e_MAF <- estimate.MAF(genodata, Nind = Nind, rho = rho, check = FALSE)
  e_b <- estimate.b(genodata, e.MAF = tmpe_e_MAF, Nind = Nind)
  if(is.null(.RM_Mat)) .RM_Mat <- haplo2RM(Allhaplo)
  ## haplo2RM can be memoised
  restmp0 <- spams.lasso(X = e_b, D = .RM_Mat$phi, return_reg_path = FALSE,
                         lambda1 = 0 , mode = 'L2ERROR', pos = TRUE)
  tmpresult <- restmp0/.RM_Mat$tau
  p <- as.numeric(tmpresult/sum(tmpresult))
  return(p)
  # if (is.null(htmat)) return(freq2haplofreq(p)) else 
  # return(as.haplofreq(htmat, p))
}


CSHAPEM <- function(genodata, rho = 0, HWE = TRUE, Allhaplo = NULL, truncthres2 = 1e-3, inferphase = FALSE,
                       udata = NULL, .RM_Mat = NULL){
  # StartTime <- Sys.time()
  q <- ncol(genodata)
  if (is.null(Allhaplo)) Allhaplo <- haplo(q)
  p_initial <- CSHAP0(genodata = genodata, Nind = 1, 
                      rho = rho, HWE = HWE, 
                      Allhaplo = Allhaplo, .RM_Mat = .RM_Mat)
  # gc(reset = TRUE)
  EM_out <- fEM(genodata = genodata, p_initial = p_initial,
                   Allhaplo = Allhaplo, truncthres2 = truncthres2, inferphase = inferphase, udata = udata)
  
  return(EM_out)
}

