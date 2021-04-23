# Sourceable

### 用来存储Phase 的类
# 
# 
# 是一个长度为T的列表，第i个元素代表了第i个个体的相型不确定性分布（是一个 data.frame ）
# 有一个 q 属性
# df结构：h1, h2, Prob。其中h1 h2 是 htc/SNPc

# Phase -------------------------------------------------------------------


Phase <- function(X, q = NULL){
  if (is.null(q)) q <- str_length(X[[1]][1,1])
  return(structure(X, class = 'Phase', q = q))
}


as.Phase <- function(...){
  UseMethod('as.Phase')
}

as.Phase.HtPhase <- function(X_HtP){
  if (is.null(X_HtP$ht$htc)) htc <- as.SNPc.SNPv(X_HtP$ht$ht) else htc <- X_HtP$ht$htc
  
  as.P.HtP.sub <- function(o){
    o2 <- matrix(htc[o[,1:2]], ncol = 2)
    colnames(o2) <- c('h1', 'h2')
    return(data.frame(o2, Prob = o[,3], stringsAsFactors = FALSE))
  }
  
  X <- llply(.data = X_HtP$Phase, .fun = as.P.HtP.sub)
  q <- X_HtP$ht$q
  return(Phase(X = X, q = q))
}

as.Phase.SpPhase <- function(X_SP, q = NULL){
  if (is.null(q)) q <- attributes(X_SP)$q
  Allhaplo <- haplo(q)
  
  as.P.SP.sub <- function(o){
    o2 <- data.frame(
      stringsAsFactors = FALSE,
      h1 = unclass(as.SNPc.SNPv(Allhaplo[, o[, 1], drop = FALSE])),
      h2 = unclass(as.SNPc.SNPv(Allhaplo[, o[, 2], drop = FALSE])),
      Prob = o[, 3]
    )
    return(o2)
  }
  X <- llply(.data = X_SP, .fun = as.P.SP.sub)
  return(Phase(X = X, q = q))
}




# HtPhase -----------------------------------------------------------------


# 
# 类似 Phase，但是多了一个 ht 组件代表 Haploset。并且相型信息存在 Phase内，此时 df 中的h1 h2 改为数字。
HtPhase <- function(ht, Phase){
  return(structure(list(ht = ht, Phase = Phase), class = 'HtPhase'))
}

as.HtPhase <- function(...){
  UseMethod('as.HtPhase')
}

as.HtPhase.Phase <- function(X_P, q = NULL){
  if (is.null(q)) q <- str_length(X[[1]][1,1])
  htc <- sort(unique(unlist(llply(.data = X_P, .fun = function(o){
    o[,1:2]
  }))))
  ht <- haploset(q = q, htc = htc)
  
  as.HtP.P.sub <- function(o){
    o2 <- matrix(match(as.matrix((o[,1:2])), htc), ncol = 2)
    colnames(o2) <- c('h1', 'h2')
    return(cbind(o2, Prob = o[,3]))
  }
  
  phase <- llply(.data = X_P, .fun = as.HtP.P.sub)
  return(HtPhase(ht = ht, Phase = phase))
}

as.HtPhase.SpPhase <- function(X_SP, q = NULL){
  if (is.null(q)) q <- attributes(X_SP)$q
  htn_all <- sort(unique(unlist(llply(
    .data = X_SP,
    .fun = function(o) {
      unlist(o[, 1:2])
    }
  ))))
  
  ht <- as.haploset(as.SNPv.integer(htn_all - 1, q = q))
  as.HtP.SP.sub <- function(o){
    o2 <- cbind(h1 = match(o[, 1], htn_all), h2 = match(o[, 2], htn_all), Prob = o[, 3])
    return(o2)
  }
  Phase <- llply(.data = X_SP, .fun = as.HtP.SP.sub)
  return(HtPhase(ht = ht, Phase = Phase))
}

as.HtPhase.SpPhase2 <- function(X_SP, Sp_Sidx, ht, q){
  ## 为了效率考虑，这里假设知道S_idx，同时强制要求q
  ## 假设 Allhaplo[,Sp_Sidx] == ht$ht
  Q_idx <- integer(2 ^ q)
  Q_idx[Sp_Sidx] <- seq_along(Sp_Sidx)
  PMat <- llply(.data = X_SP, .fun = function(o){
    o[,1] <- Q_idx[o[,1]]
    o[,2] <- Q_idx[o[,2]]
    return(o)
  })  
  return(HtPhase(ht = ht, Phase = PMat))
}

simplify.HtPhase <- function(X_HtP){
  R_idx <- rank(X_HtP$ht$htc)
  X_NewHf <- simplify.haplofreq(X_HtP$ht)
  X_NewHtP_Phase <- 
   llply(X_HtP$Phase, function(o){
   o2 <-  cbind(h1 = R_idx[o[,1]], h2 = R_idx[o[,2]], Prob = o[,3])
   return(o2)
  })
  return(HtPhase(ht = X_NewHf, Phase = X_NewHtP_Phase))
}
# SpPhase -----------------------------------------------------------------

# SpPhase
# 
# 多了一个 q 属性，此时 df 中的h1 h2 是htn ，也就是 SNPn+1。
SpPhase <- function(X, q){
 return(structure(X, class = 'SpPhase', q = q))
}

as.SpPhase <- function(...){
  UseMethod('as.SpPhase')
}

as.SpPhase.HtPhase <- function(X_HtP){
  htnlist <- as.htn.SNPv(X_HtP$ht$ht)
  q <- X_HtP$ht$q
  as.Sp.HtP.sub <- function(o){
    o2 <- matrix(htnlist[as.matrix(o[,1:2])], ncol = 2)
    colnames(o2) <- c('h1', 'h2')
    return(cbind(o2, Prob = o[,3]))
  }
  X <- llply(.data = X_HtP$Phase, .fun = as.Sp.HtP.sub)
  return(SpPhase(X = X, q = q))
}


as.SpPhase.Phase <- function(X_P, q = NULL){
  if (is.null(q)) q <- str_length(X[[1]][1,1])
  htc_all <- llply(
    .data = X_P,
    .fun = function(o) {
     o[, 1:2]
    }
  )
  htc <- sort(unique(unlist(htc_all)))
  htnlist <- as.htn.SNPc(htc)
  
  as.SP.P.sub <- function(o){
    o2 <-
      cbind(h1 = htnlist[match(o[, 1], htc)],
            h2 = htnlist[match(o[, 2], htc)],
            Prob = o[, 3])
    return(o2)
  }
  X <- llply(.data = X_P, .fun = as.SP.P.sub)
  return(SpPhase(X = X, q = q))
}


# M_Step ------------------------------------------------------------------

M.Step.SpPhase <- function(X_SP, q = NULL){
  if (is.null(q)) q <- attributes(X_SP)$q
  
}
