#include <Rcpp.h>
// #include <RcppArmadillo.h>
//#include <iostream>
#include <math.h>


using namespace Rcpp;
// using namespace RcppArmadillo;
// using namespace arma;
using namespace std;


// [[Rcpp::export]]
// 将一个整数转化为二进制形式
IntegerVector as_SNPv(int x, int q) {
  IntegerVector out(q);
  for (size_t i = 0; x > 0 & i< q; i++){
    out[q - 1 - i] = x % 2;
    x = x >> 1;
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector as_LSNPv(long int x, int q) {
  IntegerVector out(q);
  for (int i = 0; x > 0 & i< q; i++){
    out[q - 1 - i] = x % 2;
    x = x >> 1;
  }
  return out;
}

// haplo()
// [[Rcpp::export]]
IntegerMatrix haplo(int q){
  int r = pow(2, q);
  IntegerMatrix res(q, r);
  int x;
  for (int i = 0; i < r ;i++){
    x = i;
    for(int j = 0; x > 0 & j<q; j++){
      res(q-1-j,i) = x%2;
      x = x>>1;
    }
  }
  return res;
}

// 给定一系列htn，求出对应的haploset矩阵
// [[Rcpp::export]]
IntegerMatrix haplo2(int q, IntegerVector htn){
  int m = htn.length();
  IntegerMatrix res(q, m);
  int x;
  for (int i = 0; i < m; i++){
    x = htn[i] - 1 ;
    for(int j = 0; x > 0 & j<q; j++){
      res(q-1-j,i) = x%2;
      x = x>>1;
    }
  }
  return res;
}


// psi = psi1 + haplo + psi3
IntegerMatrix create_psi3(IntegerMatrix hap){
  int q = hap.nrow();
  IntegerMatrix psi3(q*(q-1)/2, hap.ncol());
  int crow = 0;
  for (int i = 1 ; i <= q - 1 ;i++){
    for (int j = i + 1; j <= q ; j++){
      for (int k=0;k<hap.ncol();k++) psi3(crow,k) = (hap(i-1,k)==0 || hap(j-1,k)==0?0:1);
      ++crow ;
    }
  }
  return psi3;
}

// psi
// 从haplo(q) 矩阵直接生成psi，节约时间&方便扩展
IntegerMatrix haplo2psi(IntegerMatrix hap){
  int q = hap.nrow();
  int r = hap.nrcol();
  IntegerMatrix psi3 = create_psi3(hap);
  IntegerMatrix result(1 + q + psi3.nrow(), hap.ncol());
  for(int i =0; i<r;i++) result(0, i) = 1;
  for (int i = 1; i < 1 + q + psi3.nrow(); i ++){
    if (i < 1 + q){
      result(i, _) = hap(i - 1, _);
      }else{
      result(i, _) = psi3(i - 1 - q, _);
    }
  }
  return result;
}

IntegerMatrix create_psi(int q){
  int r = pow(2,q);
  int row2 = q*(q-1)/2;
  int x;
  int crow=1+q;
  IntegerMatrix result(1 + q + row2, r);
  // 1
  for(int i = 0; i<r;i++) result(0, i) = 1;
  // H
  for (int i = 0; i < r ;i++){
    x = i;
    for(int j = 0; x > 0 & j<q; j++){
      result(q-j,i) = x%2;
      x = x>>1;
    }
  }
  // H^H
  for (int i = 1 ; i <= q  ;i++){
    for (int j = i + 1; j <= q ; j++){
      for (int k=0;k<r;k++) result(crow,k) = (result(i,k)==0 || result(j,k)==0?0:1);
      ++crow ;
    }
  }
  return result;
}

IntegerMatrix create_psi_slow(int q){
  IntegerMatrix hap = haplo(q);
  IntegerMatrix result = haplo2psi(hap);
  return result;
}


// 计算一个向量x的2范数
// [[Rcpp::export]]
double l2norm(NumericVector x){
  NumericVector x2 = x * x;
  double out = sum(x2);
  out = sqrt(out);
  return out;
}

// 计算一个矩阵A每一列的2范数
// [[Rcpp::export]]
NumericVector l2normM(NumericMatrix A){
  int Acol = A.ncol();
  NumericVector out(Acol);
  for (int i = 0; i < Acol; i ++){
    out[i] = l2norm(A(_, i ));
  }
  return out;
}


double l2norm_i(IntegerVector x){
  IntegerVector x2 = x * x;
  double out = sum(x2);
  out = sqrt(out);
  return out;
}

NumericVector l2normM_i(IntegerMatrix A){
  int Acol = A.ncol();
  NumericVector out(Acol);
  for (int i = 0; i < Acol; i ++){
    out[i] = l2norm_i(A(_, i ));
  }
  return out;
}

// create_RM()
// [[Rcpp::export]]
List create_RM(int q){
  IntegerMatrix H = haplo(q);
  IntegerMatrix psi = haplo2psi(H);
  NumericVector tau = l2normM_i(psi);
  int psicol = psi.ncol(), psirow = psi.nrow();
  NumericMatrix phi(psirow, psicol);
  IntegerVector psii(psirow);
  for (int i = 0;i < psicol; i++){
    psii = psi(_, i);
    phi(_, i) = as<NumericVector>(psii) / tau[i];
  }
  List result = List::create(Named("q") = q,
                             // Named("H") = H,
                             // Named("psi") = psi,
                             Named("tau") = tau,
                             Named("phi") = phi);
  return result;

}

// 允许从任意单体型配置矩阵生成RM阵
// [[Rcpp::export]]
List haplo2RM(IntegerMatrix H){
  IntegerMatrix psi = haplo2psi(H);
  NumericVector tau = l2normM_i(psi);
  int psicol = psi.ncol(), psirow = psi.nrow();
  NumericMatrix phi(psirow, psicol);
  IntegerVector psii(psirow);
  for (int i = 0;i < psicol; i++){
    psii = psi(_, i);
    phi(_, i) = as<NumericVector>(psii) / tau[i];
  }
  List result = List::create(Named("q") = H.nrow(),
                          // Named("H") = H,
                          // Named("psi") = psi,
                             Named("tau") = tau,
                             Named("phi") = phi);
  return result;

}

//////////// 优化 uniquedata

// 计算向量ux中的每个元素在x中出现的个数
// [[Rcpp::export]]
IntegerVector countint(IntegerVector ux, IntegerVector x){
  int uxl = ux.size(), xl = x.size();
  IntegerVector res(uxl);
  LogicalVector matchres(xl);
  for (int i = 0; i < ux.size(); i++){
    matchres = ux[i] == x;
    res[i] = sum(matchres);
  }
  return res;
}

// 计算向量ux中的每个元素在x中出现的个数，以及ux的每一个元素在x中哪些位置出现
// [[Rcpp::export]]
List countint_detail(IntegerVector ux, IntegerVector x){
  int uxl = ux.size(), xl = x.size();
  IntegerVector res(uxl);
  LogicalVector matchres(xl);
  IntegerVector matchidx = seq_along(x);
  List res2(ux.size());
  List result;
  for (int i = 0; i < ux.size(); i++){
    matchres = ux[i] == x;
    res2[i] = matchidx[matchres];
    res[i] = sum(matchres);
  }
  result = List::create(_["count"] = res, _["matchidx"] = res2) ;
  return result;
}


// [[Rcpp::export]]
// 判断向量x与矩阵A所有行是否相等，返回的是一个逻辑向量
LogicalVector match_row(IntegerVector x, IntegerMatrix A){
    LogicalVector res(A.nrow());
      for (int i = 0; i < A.nrow(); i++){
        res[i] = is_true(all(x == A(i,_)));
      }
   return res;
}


// [[Rcpp::export]]
// 向量x与矩阵A所有行有多少行相等，返回的是一个数
int count_row(IntegerVector x, IntegerMatrix A){
  int out;
  LogicalVector matchres = match_row(x, A);
  out = sum(matchres);
  return out;
}

// [[Rcpp::export]]
// 矩阵uA的每一行与矩阵A有多少个相等的。
IntegerVector count_row_mat(IntegerMatrix uA, IntegerMatrix A){
  IntegerVector out(uA.nrow());
  for (int i = 0; i < uA.nrow(); i++){
    LogicalVector matchres = match_row(uA(i,_), A);
    out[i] = sum(matchres);
  }
  return out;
}


// [[Rcpp::export]]
// 判断向量x与矩阵A中的哪一行相等，注：必须保证有且仅有唯一一行相等，返回的是r风格的行号坐标
int which_row(IntegerVector x, IntegerMatrix A){
  int i = 0;
  while(is_false(all(x == A(i,_)))){
    i++;
  }
  return i+1;
}

// [[Rcpp::export]]
// 矩阵uA的每一行与矩阵A有多少个相等的。并且这些相等行在A中的行坐标
List count_row_mat_detail(IntegerMatrix uA, IntegerMatrix A){
  IntegerVector out(uA.nrow());
  IntegerVector matchidx = seq_len(A.nrow());
  List res2(uA.nrow());
  List result;
  for (int i = 0; i < uA.nrow(); i++){
    LogicalVector matchres = match_row(uA(i,_), A);
    out[i] = sum(matchres);
    res2[i] = matchidx[matchres];
  }
  result = List::create(_["count"] = out, _["matchidx"] = res2) ;
  return result;
}


// [[Rcpp::export]]
// 判断某个向量x与矩阵A所有列是否相等，返回的是一个逻辑向量
LogicalVector match_col(IntegerVector x, IntegerMatrix A){
  LogicalVector res(A.ncol());
  for (int i = 0; i < A.ncol(); i++){
    res[i] = is_true(all(x == A(_,i)));
  }
  return res;
}

// [[Rcpp::export]]
// 向量x与矩阵A所有行有多少列相等，返回的是一个数
int count_col(IntegerVector x, IntegerMatrix A){
  int out;
  LogicalVector matchres = match_col(x, A);
  out = sum(matchres);
  return out;
}

// [[Rcpp::export]]
// 矩阵uA的每一列与矩阵A有多少个相等的。
IntegerVector count_col_mat(IntegerMatrix uA, IntegerMatrix A){
  IntegerVector out(uA.ncol());
  for (int i = 0; i < uA.ncol(); i++){
    LogicalVector matchres = match_col(uA(_,i), A);
    out[i] = sum(matchres);
  }
  return out;
}

// [[Rcpp::export]]
// 判断向量x与矩阵A中的哪一列相等，注：必须保证有且仅有唯一一列相等，返回的是r风格的列号坐标
int which_col(IntegerVector x, IntegerMatrix A){
  int i = 0;
  while(is_false(all(x == A(_,i)))){
    i++;
  }
  return i+1;
}

// [[Rcpp::export]]
// 矩阵uA的每一列与矩阵A有多少个相等的。并且这些相等列在A中的列坐标
List count_col_mat_detail(IntegerMatrix uA, IntegerMatrix A){
  IntegerVector out(uA.ncol());
  IntegerVector matchidx = seq_len(A.ncol());
  List res2(uA.ncol());
  List result;
  for (int i = 0; i < uA.ncol(); i++){
    LogicalVector matchres = match_col(uA(_,i), A);
    out[i] = sum(matchres);
    res2[i] = matchidx[matchres];
  }
  result = List::create(_["count"] = out, _["matchidx"] = res2) ;
  return result;
}


// [[Rcpp::export]]
// # 一个向量a与向量b连接
IntegerVector ligation_vec(IntegerVector a, IntegerVector b){
  int la = a.size(), lb = b.size();
  IntegerVector out(la + lb);
  for (int i = 0; i < la; i++){
    out(i) = a(i);
  }
  for (int i = 0; i < lb; i++){
    out(i + la) = b(i);
  }
  return out;
}


//[[Rcpp::export]]
//##  将一个矩阵重复n次，每一列分别重复n次，再并起来
IntegerMatrix rep1_mat(IntegerMatrix A, int n){
  int nA = A.nrow(), cA = A.ncol();
  IntegerMatrix out(nA, cA*n);
  for (int i = 0; i < cA; i++){
    for (int j = 0; j < n; j++){
      out(_, i * n + j) = A(_,i);
    }
  }
  return out;
}

//[[Rcpp::export]]
// Matrix is a vector with dim recording its ncol and nrow
IntegerMatrix rep1_mat_fast(IntegerMatrix A, int n){
  int rA = A.nrow(), cA = A.ncol();
  int sA = rA*cA;
  IntegerMatrix out(rA, cA*n);
  for(size_t i=0; i< rA*cA*n; i++){
    out[i] = A[i%sA];
  }
  return out;
}

