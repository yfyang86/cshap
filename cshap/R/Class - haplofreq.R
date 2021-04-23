# haplofreq0 Class ---------------------------------------------------------
# 
# 用于存储q较大时候的单体型频率，简称：hf
# 
# 必备组件：
# 
# - q：位点个数
# - ht：类似SNPv，但是带有列名，方便识别
# - p：在支撑集 ht 上的频率
# 
haplofreq0 <- function(q = NULL, ht, p, ht.rename = FALSE, check = FALSE){
  if (is.null(q)) q <- nrow(ht)
  if (ht.rename){
    colnames(ht) <- 1:ncol(ht)
  }
  if (check){
    idx <- p > 0
    p <- p[idx]
    ht <- ht[, idx, drop = FALSE]
  }
  out <- structure(list(q = q, ht = ht, p = p), class = c('haplofreq0','haploset0'))
  return(out)
}

is.haplofreq0 <- function(X){
  return('haplofreq0' %in% class(X))
}


print.haplofreq0 <- function(X.hf){
  print.haplofreq(as.haplofreq.haplofreq0(X.hf))
}


as.haplofreq0 <- function(...){
  UseMethod('as.haplofreq0')
}

as.haplofreq0.Sparsep <- function(X.sp){
  q <- X.sp$q
  ht <- t(as.matrix(X.sp$s.summary[,1:q]))
  p <- X.sp$spar_p
  return(haplofreq0(q = q, ht = ht, p = p))
}

as.haplofreq0.numeric <- function(p){
  r <- length(p)
  q <- log(r, base = 2)
  su_idx <- which(p > 0)
  ht <- as.SNPv.integer(Xn = su_idx - 1, q = q)
  return(haplofreq0(q = q, ht = ht, p = p[su_idx]))
}

# haplofreq Class ---------------------------------------------------------
# 
# 在 haplofreq0 的基础上，多了 htc 组件。
# 
# 

haplofreq <- function(q = NULL, ht = NULL, htc = NULL, p, ht.rename = FALSE, check = FALSE){
  if (is.null(q)) q <- nrow(ht)
  if (is.null(ht) && !is.null(htc)) ht <- as.SNPv.SNPc(htc)
  if (ht.rename){
    colnames(ht) <- 1:ncol(ht)
  }
  if (check){
    idx <- p > 0
    p <- p[idx]
    ht <- ht[, idx, drop = FALSE]
  }
  if (is.null(htc) && !is.null(ht)) htc <- apply(ht, 2, str_c, collapse = '')
  out <- structure(list(q = q, ht = ht, htc = htc, p = p), class = c('haplofreq','haploset'))
  return(out)
}

is.haplofreq <- function(X){
  return('haplofreq' %in% class(X))
}

print.haplofreq <- function(X.hf){
  tmp <- tbl_df(data.frame(stringsAsFactors = FALSE, Haplotype = X.hf$htc, freq = X.hf$p))
  rownames(tmp) <- NULL
  print(tmp)
}


as.haplofreq <- function(...){
  UseMethod('as.haplofreq')
}

as.haplofreq.haplofreq0 <- function(X.hf){
  X.hf$htc <- apply(X.hf$ht, 2, str_c, collapse = '')
  class(X.hf) <- c('haplofreq','haploset')
  return(X.hf)
}

as.haplofreq.Sparsep <- function(X.sp){
  q <- X.sp$q
  ht <- t(as.matrix(X.sp$s.summary[,1:q]))
  p <- X.sp$spar_p
  htc <- apply(ht, 2, str_c, collapse = '')
  return(haplofreq(q = q, ht = ht, htc = htc, p = p))
}

as.haplofreq.numeric <- function(p){
  return(as.haplofreq.haplofreq0(as.haplofreq0.numeric(p)))
}

as.numeric.haplofreq0 <- function(X.hf){
  q <- X.hf$q
  if (q > 29) return(NULL)
  base_of_2 <- 2 ^ ((q - 1):0)
  p <- numeric(2 ^ q)
  idX.hf <- c(base_of_2 %*% X.hf$ht) + 1
  p[idX.hf] <- X.hf$p
  return(p)
}

as.numeric.haplofreq <- as.numeric.haplofreq0 

as.sparseVector.haplofreq0 <- function(X.hf){
  q <- X.hf$q
  if (q > 59) return(NULL)
  base_of_2 <- 2^((q-1):0)
  idX.hf <- c(base_of_2 %*% X.hf$ht) + 1
  sp <- sparseVector(x = X.hf$p, i = idX.hf, length = 2^q) 
  return(sp)
}

as.sparseVector.haplofreq <- as.sparseVector.haplofreq0

as.haplofreq.HtDiplo <- function(X_HtD){
  q <- X_HtD$ht$q
  Inds <- nrow(X_HtD$Diplo)
  x_tab <- table(X_HtD$Diplo)
  su_idx <- as.integer(names(x_tab))
  p <- as.integer(x_tab) / Inds / 2
  res <- haplofreq(q = q, ht = X_HtD$ht$ht[, su_idx, drop = FALSE], htc = X_HtD$ht$htc[su_idx], p = p, ht.rename = TRUE)
  return(res)
}

as.haplofreq.Diplo <- function(X_D){
 tmp <-  as.HtDiplo.Diplo(X_D = X_D, sort = TRUE)
 return(as.haplofreq.HtDiplo(tmp))  
}


as.numeric.HtDiplo <- function(X_HtD){
  out1 <- as.haplofreq.HtDiplo(X_HtD)
  return(as.numeric.haplofreq(out1))
}

as.haplofreq.HtPhase <- function(X_HtP){
  q <- X_HtP$ht$q
  S <- ncol(X_HtP$ht$ht)
  as.hf.HtP.sub <- function(o){
    p <- numeric(S)
    if (nrow(o) == 1){ ## 仅有一种配置
      if (o[1,1] == o[1,2]) {
        p[o[1,1]] <- 1 ## 纯合子
      } else {
        p[o[1,1:2]] <- 0.5
      }
    } else {
      Prob <- o[,3]
      p[o[, 1]] <- Prob / 2
      p[o[, 2]] <- Prob / 2
    }
    return(p)
  }
  p <- rowMeans(sapply(X_HtP$Phase, as.hf.HtP.sub))
  res <- haplofreq(q = q, ht = X_HtP$ht$ht, htc = X_HtP$ht$htc, p = p, ht.rename = TRUE, check = TRUE)  
  return(res)
}

as.numeric.HtPhase <- function(X_HtP){
  out1 <- as.haplofreq.HtPhase(X_HtP)
  return(as.numeric.haplofreq(out1))
}

simplify <- function(...) {
  UseMethod('simplify')
}


simplify.haplofreq <- function(X.hf){
  ## 删去频率为0的，并且排序
  if (!is.null(X.hf$p)) {
  S_idx <- X.hf$p > 0
  X_ht <- X.hf$ht[, S_idx, drop = FALSE]
  X_htc <- X.hf$htc[S_idx]
  X_p <- X.hf$p[S_idx]
  } else {
    X_ht <- X.hf$ht
    X_htc <- X.hf$htc
  }
  O_idx <- order(X_htc)
  out <- haplofreq(q = X.hf$q, ht = X_ht[, O_idx, drop = FALSE], htc = X_htc[O_idx], p = NULL, ht.rename = TRUE)
  if (!is.null(X.hf$p)) {
    out$p = X_p[O_idx]
  }
  return(out)
}