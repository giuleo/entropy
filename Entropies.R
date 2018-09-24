# =================================================================
# Wrapper --
# Simplified and modified crqa function (by Coco & Dale, 2014).
# It takes as an input the categorical time series, runs
# a categorical recurrence analysis and returns the various estimates 
# of entropy: dlEntr (based on diagonal lines), baEntr (based 
# on blocks area), vlEntr (based on vertical lines) and hlEntr
# (based on horizontal lines). It also returns the rEntr measure
# from the original function (normalized entropy based on diagonals) 
# =================================================================

Entropies <- function (
  ts1,
  ts2,
  radius = 0.1,
  mindiagline = 2,
  minvertline = 2,
  recpt = F)
{
  if (recpt == FALSE) {  # ts1 and ts2 are time series
    ts1 = as.vector(as.matrix(ts1))
    ts2 = as.vector(as.matrix(ts2))
    
    v1l = length(ts1)
    v2l = length(ts2)
    
    # Distance Matrix
    dm = rdist(ts1, ts2)
    
    # From Distance Matrix to Recurrence Plot
    # r and c are the row and column indices for recurrence points
    ind = which(dm <= radius, arr.ind = TRUE)
    r = ind[, 1]
    c = ind[, 2]
  }
  
  # here we enter if recpt == TRUE. 
  # ts1 in this case is a Recurrence Plot
  else { 
    ts1 = matrix(as.logical(ts1), ncol = ncol(ts1))
    v1l = nrow(ts1)
    v2l = ncol(ts1)
    ind = which(ts1 > 0, arr.ind = TRUE)
    r = ind[, 1]
    c = ind[, 2]
  }
  
  # Recurrence Quantification -
  # we enter here if there is at least one recurrence point
  if (length(r) != 0 & length(c) != 0) {
    
    S = sparseMatrix(r, c, dims = c(v1l, v2l))
    S = t(S)                      # transpose to put ts1 in columns (?)
    S = theiler(S, 0)             # theiler function in crqa package
    
    # MY PIECE OF CODE FOR CATEGORICAL ENTROPY
    
    baEntropy <- catEntropy(S)$Entropy
    vhentropy <- vhEntropy(S, minvertline)
    vlEntr <- vhentropy[[1]]
    hlEntr <- vhentropy[[2]]
    
    # extract non-zero diagonal lines in the form of a matrix B
    spdiagonalize = spdiags(S)   
    B = spdiagonalize$B
    
    # add a row of FALSEs before and after B
    if (is.vector(B)) {
      false = rep(FALSE, length(B))
      B = rbind(false, B, false, deparse.level = 0)
    }
    else {
      false = rep(FALSE, ncol(B))
      B = as.matrix(B)
      B = rbind(false, B, false, deparse.level = 0)
    }
    
    # Find and sort all the diagonal lines
    diaglines = sort(diff(which(B == FALSE)) - 1, decreasing = TRUE)
    diaglines = diaglines[-which(diaglines < mindiagline)]
    
    # If indeed there are diagonal lines compute Recurrence Measures
    if (length(diaglines) != 0) {
      tabled = as.data.frame(table(diaglines))
      total = sum(tabled$Freq)
      p = tabled$Freq / total
      del = which(p == 0)
      if (length(del) > 0) {
        p = p[-del]
      }
      entropy = -sum(p * log(p))
      relEntropy = entropy / (-1 * log(1 / nrow(tabled)))
    }
    else {
      entropy = NA
      relEntropy = NA
    }
    results = list(
      dlENTR = entropy,
      rENTR = relEntropy,
      baENTR = baEntropy,
      vlENTR = vlEntr,
      hlENTR = hlEntr
    )
  }
  
  else { # In the case of not a single recurrence point
    results = list(
      dlENTR = NA,
      rENTR = NA,
      baENTR = NA,
      vlENTR = NA,
      hlENTR = NA
    )
  }
  
  # What it spits out
  return(results)
}


# ---------------------------
# Computation of entropy from vector 
# ---------------------------

Shannon <- function(x){
  tab <- as.data.frame(table(x))
  par <- tab$Freq / sum(tab$Freq)
  H = -sum(par * log(par))
  return(H)
}
