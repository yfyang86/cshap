# 目标是定义三种双体型类：
#  Diplo 应该尽量原始，与genodata矩阵对应，使用两个行向量代表双体型，适用于任何规模的q和T。
#  Diploc 是Diplo的字符串版本，Tx2矩阵
#  SparseDp应该与 Sparsep 对应，使用两个（排序的）数字代表个体的双体型，适用于q较小的时；
#  HtDiplo 应该与 haplofreq 对应，使用 haploset 的数字代表个体双体型，适用于q中等规模并且T较大时；


# Diplo Class -------------------------------------------------------------
# 2T x q 矩阵，每一行为一个单体型向量，第 2i-1 和 2i 行代表第 i 个个体的双体型配置
Diplo <- function(DpMat){
  Inds <- nrow(DpMat)/2
  oldrowname <- str_c('IND', 1:Inds)
  rownames(DpMat) <- c(rbind(str_c(oldrowname, '_1'), str_c(oldrowname, '_2')))
  class(DpMat) <- 'Diplo'
  return(DpMat)
}


as.Diplo <- function(...){
  UseMethod('as.Diplo')
}

# SpDiplo 2 Diplo
as.Diplo.SpDiplo <- function(X_SpD, Allhaplo = NULL, q = NULL){
  if (is.null(q)) q <- attributes(X_SpD)$q
  if (is.null(Allhaplo)) Allhaplo <- haplo(q)
  Inds <- nrow(X_SpD)
  htns <- c(t(X_SpD))
  res <- t(Allhaplo[, htns, drop = FALSE])
  if (!is.null(rownames(X_SpD))) {
    oldrowname <- rownames(X_SpD)
  } else {
    oldrowname <- str_c('IND', 1:Inds)
  }
  rownames(res) <- c(rbind(str_c(oldrowname, '_1'), str_c(oldrowname, '_2')))
  class(res) <- 'Diplo'
  return(res)
}

as.Diplo.HtDiplo <- function(X_HtD){
  idx <- c(t(X_HtD$Diplo))
  res <- t(X_HtD$ht$ht[,idx, drop = FALSE])
  return(Diplo(res))
}



as.Diplo.Diploc <- function(X_Dc){
  q <- str_length(X_Dc[1,1])
  Inds <- nrow(X_Dc)
  htc_all <- unique(X_Dc)
  htmat <- as.SNPvt.SNPc(htc_all)
  res <- matrix(NA, nrow = 2*Inds, ncol = q)
  
  res[(1:Inds)*2 - 1,] <- htmat[match(X_Dc[, 1], htc_all), ]
  res[(1:Inds)*2,] <- htmat[match(X_Dc[, 2], htc_all), ]
  return(Diplo(res))
}
# Phase 2 Diplo
as.Diplo.Phase <- function(X_P){
  as.D.P.sub <- function(o){
    as.SNPvt.SNPc(o[which.max(o[,3]), 1:2])
  }
  
  DpMat <- do.call(rbind, llply(.data = X_P, .fun = as.D.P.sub))
  return(Diplo(DpMat = DpMat))
}


# Diploc Class ------------------------------------------------------------
Diploc <- function(X){
  class(X) <- 'Diploc'
  return(X)
}

as.Diploc <- function(...){
  UseMethod('as.Diploc')
}

as.Diploc.Diplo <- function(X_D) {
  Inds <- nrow(X_D) / 2
  res <- cbind(
    h1 = apply(
      X_D[(1:Inds) * 2 - 1,],
      MARGIN = 1,
      FUN = str_c,
      collapse = ''
    ) ,
    h2 = apply(
      X_D[(1:Inds) * 2,],
      MARGIN = 1,
      FUN = str_c,
      collapse = ''
    )
  )
  rownames(res) <- str_c('IND', 1:Inds)
  class(res) <- 'Diploc'
  return(res)
}  


as.Diploc.SpDiplo <- function(X_SpD, q = NULL, Allhaplo = NULL) {
  if (is.null(q)) q <- attributes(X_SpD)$q
  if (is.null(Allhaplo)) {
    Allhaplo <- haplo(q)
  }
  htn_all <- unique(c(X_SpD))
  htc <- as.SNPc.SNPv(Allhaplo[, htn_all, drop = FALSE])
  
  res <- cbind(
    htc[match(X_SpD[,1], htn_all)],
    htc[match(X_SpD[,2], htn_all)]
  )
  
  class(res) <- 'Diploc'
  return(res)
}


as.Diploc.HtDiplo <- function(X_HtD) {
  htc <- X_HtD$ht$htc
  if (is.null(htc)) {
    htc <- as.SNPc.SNPv(X_HtD$ht$ht)
  }
  
  res <- cbind(htc[X_HtD$Diplo[, 1]],
               htc[X_HtD$Diplo[, 2]])
  class(res) <- 'Diploc'
  return(res)
}



as.Diploc.Phase <- function(X_P, q = NULL) {
  if (is.null(q)) q <- attributes(X_P)$q
  res <- laply(
    .data = X_P,
    .fun = function(o) {
      unlist(o[which.max(o[, 3, drop = TRUE]), 1:2])
    }
  )
  return(Diploc(res))
}


# SpDiplo Class -----------------------------------------------------------
# Dpdata: Tx2矩阵，每一行代表一个个体的两个单体型在Haplo(q)中对应的htn
SpDiplo <- function(Dpdata, q){
  res <- structure(Dpdata, class = 'SpDiplo', q = q)
  return(res)
}

as.SpDiplo <- function(...){
  UseMethod('as.SpDiplo')
}

as.SpDiplo.Diplo <- function(X_D){
  q <- ncol(X_D)
  Inds <- nrow(X_D)/2
  base2 <- 2^((q:1) - 1)
  htns <-  X_D %*% base2 + 1
  X2 <- cbind(htns[2*(1:Inds) - 1], htns[2*(1:Inds)])
  rownames(X2) <- str_sub(rownames(X_D)[2*(1:Inds)], 1, -3)
  return(SpDiplo(X2, q))
}

as.SpDiplo.HtDiplo <- function(X_HtD){
  q <- X_HtD$ht$q
  base2 <- 2^((q:1) - 1)
  htns <- base2 %*% X_HtD$ht$ht + 1
  X2 <- cbind(htns[X_HtD$Diplo[,1]], htns[X_HtD$Diplo[,2]])
  rownames(X2) <- rownames(X_HtD$Diplo)
  return(SpDiplo(X2, q))
}


as.SpDiplo.Diploc <- function(X_Dc) {
  htc_all <- unique(c(X_Dc))
  htn_all <- as.htn.SNPc(htc_all)
  q <- str_length(X_Dc[1, 1])
  res <- cbind(htn_all[match(X_Dc[, 1], htc_all)],
               htn_all[match(X_Dc[, 2], htc_all)])
  return(SpDiplo(Dpdata = res, q = q))
}


as.SpDiplo.SpPhase <- function(X_SP, q = NULL) {
  if (is.null(q)) q <- attributes(X_SP)$q
  res <-  laply(
    .data = X_SP,
    .fun = function(o) {
      unlist(o[which.max(o[, 3, drop = TRUE]), 1:2])
    }
  )
  return(SpDiplo(Dpdata = res, q = q))
}


# HtDiplo Class-----------------------------------------------------------------
# ht: A haplofreq Or haploset or haploset0
# Diplo: Tx2 矩阵，每一行代表一个体
HtDiplo <- function(ht, Diplo){
  res <- structure(list(ht = ht, Diplo = Diplo), class = 'HtDiplo')
  return(res)
}

# 化简HtDiplo，删除Ht中未出现在Diplo中的单体型，并且对Ht矩阵进行排序。
simplify.HtDiplo  <-  function(X_HtD){
  R_idx <- rank(X_HtD$ht$htc)
  O_idx <- order(X_HtD$ht$htc) 
  X_HtD$ht$htc <- X_HtD$ht$htc[O_idx]
  X_HtD$ht$ht <- X_HtD$ht$ht[, O_idx, drop = FALSE]
  if (!is.null(X_HtD$ht$count)) X_HtD$ht$count <- X_HtD$ht$count[O_idx]
  if (!is.null(X_HtD$ht$p))X_HtD$ht$p <- X_HtD$ht$p[O_idx]
  New_D <- cbind(R_idx[X_HtD$Diplo[, 1]], R_idx[X_HtD$Diplo[, 2]])
  X_HtD$Diplo <- New_D
  return(X_HtD)
}

as.HtDiplo <- function(...){
  UseMethod('as.HtDiplo')
}

as.HtDiplo.SpDiplo <- function(X_SpD, Allhaplo = NULL, q = NULL){
  if (is.null(q)) q <- attributes(X_SpD)$q
  if (is.null(Allhaplo)) Allhaplo <- haplo(q)
  uhtn <- sort(unique(X_SpD))
  ht <- as.haploset0(Allhaplo[, uhtn, drop = FALSE])
  Diplo <- cbind(match(X_SpD[,1], uhtn), match(X_SpD[,2], uhtn))
  rownames(Diplo) <- rownames(X_SpD)
  return(HtDiplo(ht = ht, Diplo = Diplo))
}


as.HtDiplo.Diplo <- function(X_D, sort = FALSE){
  # sort 表示生成的ht中是否需要排序
  # 如果需要排序，先变成 Diploc 更好，因为向量无法排序
  if (sort){
    X_Dc <- as.Diploc.Diplo(X_D)
    return(as.HtDiplo.Diploc(X_Dc))
  }
  
  ug <- uniquegenodata(unclass(X_D), MARGIN = 1, detail = TRUE)
  ht <- as.haploset0(t(ug$udata))
  q <- ncol(X_D)
  Inds <- nrow(X_D)/2
  X2 <- cbind(ug$uidx[ 2*(1:Inds) - 1], ug$uidx[ 2*(1:Inds)])
  rownames(X2) <- str_sub(rownames(X_D)[2*(1:Inds)], 1, -3)
  return(HtDiplo(ht, Diplo = X2))
}


as.HtDiplo.Diploc <- function(X_Dc) {
  htc <- sort(unique(X_Dc))
  q <- str_length(X_Dc[1,1])
  ht <- haploset(q = q, htc = htc, ht.rename = TRUE)
  
  res <- cbind(match(X_Dc[, 1], htc),
               match(X_Dc[, 2], htc))
  
  return(HtDiplo(ht = ht, Diplo = res))
}

as.HtDiplo.HtPhase <- function(X_HtP) {
  res <- laply(X_HtP$Phase, function(o) {
    unlist(o[which.max(o[, 3, drop = TRUE]), 1:2])
  })
  return(HtDiplo(ht = X_HtP$ht, Diplo = res))
}

# 从单体型配置转化为基因型 ------------------------------------------------------------
as.genodata <- function(...){
  UseMethod('as.genodata')
}

as.genodata.Diplo <- function(X){
  Inds <- nrow(X)/2
  res <- X[2*(1:Inds) - 1, ] +  X[2*(1:Inds), ]
  if (!is.null(rownames(res))) rownames(res) <- str_sub(rownames(X)[2*(1:Inds)], 1, -3)
  return(res)
}

as.genodata.SpDiplo <- function(X, Allhaplo = NULL){
  q <- attributes(X)$q
  if (is.null(Allhaplo)) Allhaplo <- haplo(q)
  res <- t(Allhaplo[, X[,1], drop = FALSE] + Allhaplo[, X[,2], drop = FALSE])
  if (!is.null(rownames(res))) rownames(res) <- rownames(X)
  return(res)
}

as.genodata.HtDiplo <- function(X){
  ht <- X$ht$ht
  res <- t(ht[, X$Diplo[,1], drop = FALSE] + ht[, X$Diplo[,2], drop = FALSE])
  if (!is.null(rownames(res))) rownames(res) <- rownames(X$Diplo)
  return(res)
}

# UnitTest

