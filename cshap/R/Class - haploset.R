library(stringr)
library(pryr)
library(dplyr)

# haploset0 Class ---------------------------------------------------------
# 
# 用于存储任意长的单体型集。
# 
# 必备组件：
# 
# - q：位点个数
# - ht：类似SNPv，但是带有列名，方便识别。
#   q x m 矩阵，每一列是一个单体型，每一行代表一个位点

haploset0 <- function(q = NULL, ht, ht.rename = FALSE){
  # 默认不加列名，省时间
  if (is.null(q)) q <- nrow(ht)
  if (ht.rename){
    colnames(ht) <- 1:ncol(ht)
  }
  out <- structure(list(q = q, ht = ht), class = 'haploset0')
  return(out)
}

as.haploset0 <- function(...){
  UseMethod('as.haploset0')
}

as.haploset0.SNPv <- function(Xv){
  return(haploset0(ht = Xv))
}

as.haploset0.SNPn <- function(Xn){
  return(haploset0(ht = as.SNPv.SNPn(Xn)))
}

as.haploset0.SNPc <- function(Xc){
  return(haploset0(ht = as.SNPv.SNPc(Xc)))
}

as.haploset0.matrix <- function(X){
  return(haploset0(ht = X))
}

as.haploset0.character <- function(X){
  return(haploset0(ht = as.SNPv.SNPc(X)))
}


# haploset Class ---------------------------------------------------------
# 
# 简称ht，比 haploset0 多一个组件：
# 
# - htc：类似SNPc
# 
# 
haploset <- function(q = NULL, ht = NULL, htc = NULL, ht.rename = FALSE){
  if (is.null(q)) q <- nrow(ht)
  if (is.null(ht) && !is.null(htc)) ht <- as.SNPv.SNPc(htc)
  if (ht.rename){
    colnames(ht) <- 1:ncol(ht)
  }
  if (is.null(htc)) htc <- apply(ht, 2, str_c, collapse = '')
  out <- structure(list(q = q, ht = ht, htc = htc), class = 'haploset')
  return(out)
}


as.haploset.haploset0 <- function(X){
  X$htc <- apply(X$ht, 2, str_c, collapse = '')
  class(X) <- 'haploset'
  return(X)
}

as.haploset <- function(X){
  X0 <- as.haploset0(X)
  X0$htc <- apply(X0$ht, 2, str_c, collapse = '')
  class(X0) <- 'haploset'
  return(X0)
}


# 草稿纸 ---------------------------------------------------------------------

if (FALSE){
  haploset0(ht = ht, q = NULL, ht.name = TRUE)
  haploset(ht = ht)
}
