# Sourceable



as.sparseVector <- function(...){
  UseMethod('as.sparseVector')
}



as.sparseVector.numeric <- function(p){
  # 假设 p 是一个向量
  p_idx <- which(p > 0)
  px <- p[p_idx] 
  return(sparseVector(x = px, i = p_idx, length = length(p)))
}



# 描述单体型频率的类 Sparsep -------------------------------------------------------
# 用于存储当q较小时候的单体型频率
# 
# 必备组件：
# 
# - q：位点个数
# - p：完整的频率向量，对应单体型顺序是按照 000, 001, 101,...,111 排列的。
# - su_p：支撑集 S, htn = which(p>0) 。
# - spar_p：频率在支撑集 S 上的分量，p_S
# 
# 可选组件：
# 
# - f.summary： A tibble，记录了完整的位点信息和频率
# - s.summary：只记录了支撑集上的位点信息和频率

Sparsep <- function(q, su_p, spar_p, p = NULL, check = TRUE){
  # check：是否做检查——若为FALSE，直接包装，不保证一切
  r <- 2^q
  if (check){
    idx <- spar_p > 0
    su_p <- su_p[idx]
    spar_p <- spar_p[idx]
    oidx <- order(su_p)
    su_p <- su_p[oidx]
    spar_p <- spar_p[oidx]
  }
  if (is.null(p)){
    p <- matrix(0, nrow = r, ncol = 1)
    p[su_p] <- spar_p
  } else if (!is.matrix(p)){
    p <- matrix(p)
  }
  rownames(p) <- paste0('ht', 1:r)
  
  # 制作 Summary
  AllHaplo <- haplo(q)
  subhaplo <- t(AllHaplo[,su_p, drop = FALSE])
  colnames(subhaplo) <- paste0('loci', 1:q)
  rownames(subhaplo) <- paste0('ht', su_p)
  s.summary <- data.frame(stringsAsFactors = FALSE, subhaplo, freq = spar_p, SNPn = su_p-1, htn = su_p)
  
  out <- list(q = q, p = p, su_p = su_p, spar_p = spar_p, s.summary = s.summary)
  class(out) <- 'Sparsep'
  return(out)
}


# 方法 ----------------------------------------------------------------------
print.Sparsep <- function(X){
  print(X$s.summary)
}



# 转换到 Sparsep --------------------------------------------------------------
as.Sparsep <- function(...){
  UseMethod('as.Sparsep')
}

as.Sparsep.numeric <- function(p){
  r <- length(p)
  q <- log(r, base = 2)
  su_p <- which(p > 0)
  spar_p <- p[su_p]
  return(Sparsep(q = q, su_p = su_p, spar_p = spar_p, p = p))
}

as.Sparsep.haplofreq <- function(X.hf){
  q <- X.hf$q
  su_p <- as.SNPn.SNPv(X.hf$ht) + 1
  spar_p <- X.hf$p
  return(Sparsep(q = q, su_p = su_p, spar_p = spar_p))
}

as.Sparsep.haplofreq0 <- function(X.hf){
  q <- X.hf$q
  su_p <- as.SNPn.SNPv(X.hf$ht) + 1
  spar_p <- X.hf$p
  return(Sparsep(q = q, su_p = su_p, spar_p = spar_p))
}

as.numeric.Sparsep <- function(X.sp){
  return(c(X.sp$p))
}


as.Sparsep.SpDiplo <- function(X_SpD, q = NULL){
  if (is.null(q)) q <- attributes(X_SpD)$q
  Inds <- nrow(X_SpD)
  x_p <- table(X_SpD)
  su_p <- as.integer(names(x_p))
  spar_p <- as.integer(x_p)/(2*Inds)
  p <- numeric(2^q)
  p[su_p] <- spar_p
  res <- Sparsep(q = q, su_p = su_p, spar_p = spar_p, p = p)
  return(res)
}


as.numeric.SpDiplo <- function(X_SpD, q = NULL){
  if (is.null(q)) q <- attributes(X_SpD)$q
  Inds <- nrow(X_SpD)
  x_p <- table(X_SpD)
  su_p <- as.integer(names(x_p))
  spar_p <- as.integer(x_p)/(2*Inds)
  p <- numeric(2^q)
  p[su_p] <- spar_p
  return(p)
}

as.Sparsep.SpPhase <- function(X_SpP, q = NULL){
  if (is.null(q)) q <- attributes(X_SpP)$q
  Inds <- length(X_SpP)
  p <- rowMeans(M_Steps(E_Results = X_SpP, q = q))
  res <- as.Sparsep.numeric(p)
  return(res)
}

as.numeric.SpPhase <- function(X_SpP, q = NULL){
  if (is.null(q)) q <- attributes(X_SpP)$q
  p <- rowMeans(M_Steps(E_Results = X_SpP, q = q))
  return(p)
}

