CS_HAP <- function(genodata, N = 1, rho = 0, HWE = TRUE, phasing = FALSE){
  q <- ncol(genodata)
  if (N == 1){
    ## individual
    if (q < 20){
      ## short
      return(CSHAPEM(genodata = genodata, rho = rho, HWE = HWE, inferphase = phasing))
    } else if (q < 500){
      ## PL
      return(PLCSHAPLEM(genodata = genodata, rho = rho, HWE = HWE, inferphase = phasing))
    } else {
      ## only phasing
      return(CSHAP_2SNP(genodata = genodata, ligation_step = 1))
    }
  } else {
    ## Pooling
    if (q < 20){
      ## short
      return(CSHAPAEM(genodata = genodata, Nind = N, rho = rho, HWE = HWE))
    } else {
      return(PLCSHAPLAEM(genodata = genodata, Nind = N, rho = rho, HWE = HWE))
    }
  }

}
