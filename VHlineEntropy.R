vhEntropy <- function(X = sparseMatrix(), minvertline) 
{
  if (class(X) != "ngCMatrix") {
    stop("Input must be a Sparse Matrix")
  }
  
  vhlines <- vhline(X, minvertline)
  
  vlines <- vhlines[[1]]
  hlines <- vhlines[[2]]
  
  if (length(vlines) > 0) {
    tabvl <- as.data.frame(table(vlines))
    pvl <- tabvl$Freq / sum(tabvl$Freq)
    vlEntropy = -sum(pvl * log(pvl))
  }
  else {
    vlEntropy <- NA
  }
  
  if (length(hlines) > 0) {
    tabhl <- as.data.frame(table(hlines))
    phl <- tabhl$Freq / sum(tabhl$Freq)
    hlEntropy = -sum(phl * log(phl))
  }
  else {
    hlEntropy <- NA
  }
  
  return(list(vlEntropy = vlEntropy, hlEntropy = hlEntropy))
}



vhline <- function (x, minvertline)
{
  xv = x
  xh = t(x)
  
  xv = rBind(xv, rep(0, ncol(xv)), deparse.level = 0)
  xh = rBind(xh, rep(0, ncol(xh)), deparse.level = 0)
  
  xv = as.vector(xv)
  xh = as.vector(xh)
  
  z = diff(xv)
  w = diff(xh)
  z0 = which(z == 1)
  w0 = which(w == 1)
  z1 = which(z == -1)
  w1 = which(w == -1)
  if (z0[1] > z1[1]) {
    z0 = c(0, z0)
  }
  
  if (w0[1] > w1[1]) {
    w0 = c(0, w0)
  }
  
  if (length(z0) > length(z1)) {
    z0 = z0[-length(z0)]
  }
  
  if (length(w0) > length(w1)) {
    w0 = w0[-length(w0)]
  }
  
  t = sort(z1 - z0)
  u = sort(w1 - w0)
  t1 = t[which(t - minvertline > 0)]
  u1 = u[which(u - minvertline > 0)]
  
  return(list(vlines = t1, hlines = u1))
}