# 描述单体型的类
# 
# SNPv
# 
# 以向量（矩阵）形式存储的单体型，矩阵的每个列向量是一个单体型。
## SNPvt : 与SNPv类似，但是每一行作为一个单体型。
# 
# SNPc
# 
# 以字符串（向量）形式存储的单体型，每个字符串代表一个单体型。
# 
# SNPn
# 
# 以整数（向量）形式存储的单体型，有一个额外属性 q 存储长度，每个数字是单体型向量的字符表示。
# Sourceable
# 类的定义 --------------------------------------------------------------------
SNPn <- function(x, q){
  out <- structure(x, q = q)
  class(out) <- 'SNPn'
  return(out)
}

SNPv <- function(ht){
  class(ht) <- 'SNPv'
  return(ht)
}

SNPc <- function(xc){
  class(xc) <- 'SNPc'
  return(xc)
}

## 
as.SNPvt <- function(...){
  UseMethod('as.SNPvt')
}

as.SNPvt.SNPv <- function(Xv){
  return(unclass(t(Xv)))
}

as.SNPvt.SNPn <- function(Xn, q = NULL){
  if (is.null(q)) q <- attributes(Xn)$q
  Xv <- sapply(Xn, as_LSNPv, q = q)
  return(t(Xv))
}

as.SNPvt.SNPc <- function(Xc, q = NULL){
  Xv <- simplify2array(str_split(Xc, ''))
  storage.mode(Xv) <- 'integer'
  return(t(Xv))
}


# 相互转换 --------------------------------------------------------------------
as.SNPv <- function(...){
  UseMethod('as.SNPv')
}

as.SNPv.SNPn <- function(Xn, q = NULL){
  if (is.null(q)) q <- attributes(Xn)$q
  q <- attributes(Xn)$q
  Xv <- sapply(Xn, as_LSNPv, q = q)
  class(Xv) <- 'SNPv'
  return(Xv)
}

as.SNPv.SNPc <- function(Xc){
  Xv <- simplify2array(str_split(Xc, ''))
  storage.mode(Xv) <- 'integer'
  class(Xv) <- 'SNPv'
  return(Xv)
}

as.SNPv.htn <- function(htn, q){
  Xv <- sapply(htn - 1, as_LSNPv, q = q)
  class(Xv) <- 'SNPv'
  return(Xv)
}

as.SNPv.integer <- function(Xn, q = NULL){
  Xv <- sapply(Xn, as_LSNPv, q = q)
  class(Xv) <- 'SNPv'
  return(Xv)
}

as.SNPv.character <- function(Xc){
  Xv <- simplify2array(str_split(Xc, ''))
  storage.mode(Xv) <- 'integer'
  class(Xv) <- 'SNPv'
  return(Xv)
}

##
as.SNPn <- function(...){
  UseMethod('as.SNPn')
}

as.SNPn.SNPv <- function(Xv){
  q <- nrow(Xv)
  base2 <- 2^((q:1) - 1)
  out <- structure(c(base2 %*% Xv), q = q)
  class(out) <- 'SNPn'
  return(out)
}

as.SNPn.SNPc <- function(Xc){
  return(as.SNPn.SNPv(as.SNPv.SNPc(Xc)))
}



##
as.SNPc <- function(...){
  UseMethod('as.SNPc')
}

as.SNPc.SNPv <- function(Xv){
  res <- apply(Xv, 2, stringr::str_c, collapse = '')
  class(res) <- 'SNPc'
  return(res)
}


as.SNPc.SNPn <- function(Xn, q = NULL){
  if (is.null(q)) q <- attributes(Xn)$q
  return(as.SNPc.SNPv(as.SNPv.SNPn(Xn, q = q)))
}

as.SNPc.htn <- function(htn, q){
  return(as.SNPc.SNPv(as.SNPv.htn(htn = htn, q = q)))
}

as.SNPc.matrix <- function(Xv){
  res <- apply(Xv, 2, stringr::str_c, collapse = '')
  class(res) <- 'SNPc'
  return(res)
}

as.integer.SNPc <- function(Xc){
  Xv <- as.SNPv.SNPc(Xc)
  q <- nrow(Xv)
  base2 <- 2^((q:1) - 1)
  out <- c(base2 %*% Xv)
  return(out)
}


# htn ---------------------------------------------------------------------
as.htn <- function(...){
  UseMethod('as.htn')
}

as.htn.SNPn <- function(Xn){
  return(c(Xn) + 1)
}

as.htn.SNPv <- function(Xv){
  q <- nrow(Xv)
  if (is.null(q)) q <- length(Xv)
  base2 <- 2^((q:1) - 1)
  out <- c(base2 %*% Xv) + 1
  return(out)
}

## 区别在于，每一行是一个单体型
as.htn.SNPvt <- function(Xvt, q = NULL){
  q <- ncol(Xvt)
  if (is.null(q)) {
    q <- length(Xvt)
    base2 <- 2^((q:1) - 1)
    out <-   sum(Xvt * base2) + 1
  } else {
  base2 <- 2^((q:1) - 1)
  out <- c(Xvt %*% base2) + 1
  }
  return(out)
}


as.htn.SNPc <- function(Xc){
  Xv <- as.SNPv.SNPc(Xc)
  q <- nrow(Xv)
  base2 <- 2^((q:1) - 1)
  out <- c(base2 %*% Xv) + 1
  return(out)
}
