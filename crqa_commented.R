# =================================================================
# The original crqa function by Coco & Dale commented and modified
# I'll try to add categorical entropy
# =================================================================

crqa_gl <- function (
  ts1,
  ts2,
  delay,
  embed,
  rescale,
  radius,
  normalize = 0,
  mindiagline = 2,
  minvertline = 2,
  tw = 0,
  whiteline = F,
  recpt = F,
  side = "both",
  checkl = list(do = F))
{
  v11 = v21 = NULL
  if (recpt == FALSE) {
    ts1 = as.vector(as.matrix(ts1))
    ts2 = as.vector(as.matrix(ts2))
    
    # The following is useless it seems. ts1 and ts2 already converted in vectors
    if (is.matrix(ts1)) { 
      stop("Your data must consist of a single column of data.")
    }
    if (is.matrix(ts2)) {
      stop("Your data must consist of a single column of data.")
    }
    
    # Boh!
    if (checkl$do == TRUE) {
      tsnorm = checkts(ts2, ts1, checkl$datatype, checkl$thrshd, checkl$pad)
      if (tsnorm[[2]] == FALSE) {
        stop(
          "Time-series difference longer than threshold. Increase threshold, 
          or set checkl$do = FALSE avoiding normalization of ts"
        )
      }
      else {
        ts1 = tsnorm[[1]][, 1]
        ts2 = tsnorm[[1]][, 2]
      }
      }
    
    # More error hunting
    if (length(ts1) < embed * delay) {
      stop("Phase-space (embed*delay) longer than ts1")
    }
    if (length(ts2) < embed * delay) {
      stop("Phase-space (embed*delay) longer than ts2")
    }
    
    # Normalization
    if (normalize > 0) {
      switch(normalize, {
        1
        ts1 = (ts1 - min(ts1))
        ts1 = ts1 / max(ts1)
        ts2 = (ts2 - min(ts2))
        ts2 = ts2 / max(ts2)
      }, {
        2
        ts1 = (ts1 - mean(ts1)) / sd(ts1)
        ts2 = (ts2 - mean(ts2)) / sd(ts2)
      })
    }
    
    # Embedding (creates time delayed vector copies of the original ts)
    for (loop in 1:embed) {
      vectorstart = (loop - 1) * delay + 1
      vectorend = length(ts1) - ((embed - loop) * delay)
      assign(paste("v1", loop, sep = ""), ts1[vectorstart:vectorend])
    }
    for (loop in 1:embed) {
      vectorstart = (loop - 1) * delay + 1
      vectorend = length(ts2) - ((embed - loop) * delay)
      assign(paste("v2", loop, sep = ""), ts2[vectorstart:vectorend])
    }
    
    # Embedding (create a matrix from the time delayed copies with cbind)
    dimts1 = dimts2 = vector()
    for (loop in 1:embed) {
      if (loop == 1) {
        dimts1 = v11
      }
      else {
        eval(parse(
          text = paste(
            "dimts1 = cbind(dimts1,",
            paste("v1", loop, sep = ""),
            ", deparse.level = 0)",
            sep = ""
          )
        ))
      }
    }
    for (loop in 1:embed) {
      if (loop == 1) {
        dimts2 = v21
      }
      else {
        eval(parse(
          text = paste(
            "dimts2 = cbind(dimts2,",
            paste("v2", loop, sep = ""),
            ", deparse.level = 0)",
            sep = ""
          )
        ))
      }
    }
    
    v1l = length(v11)
    v2l = length(v21)
    
    # Distance Matrix
    dm = rdist(dimts1, dimts2)
    
    # Rescaling of the Distance matrix
    if (rescale > 0) {
      switch(rescale, {
        1
        rescaledist = mean(dm)
        dmrescale = (dm / rescaledist) * 100
      }, {
        2
        rescaledist = max(dm)
        dmrescale = (dm / rescaledist) * 100
      })
    }
    else {
      dmrescale = dm
    }
    
    # From Distance Matrix to Recurrence Plot
    # r and c are the row and column indices for recurrence points
    ind = which(dmrescale <= radius, arr.ind = TRUE)
    r = ind[, 1]
    c = ind[, 2]
}
  
  # here we enter if recpt == TRUE. 
  # But what does it mean? What does it do? Where is ts2?
  # Maybe ts1 in this case is a Recurrence Matrix already?
  else { 
    ts1 = matrix(as.logical(ts1), ncol = ncol(ts1))
    v1l = nrow(ts1)
    v2l = ncol(ts1)
    ind = which(ts1 > 0, arr.ind = TRUE)
    r = ind[, 1]
    c = ind[, 2]
  }
  
  # - Recurrence Quantification -
  # we enter here if there is at least one recurrence point
  if (length(r) != 0 & length(c) != 0) {
    
    S = sparseMatrix(r, c, dims = c(v1l, v2l))
    S = t(S)                      # transpose to put ts1 in columns (?)
    S = theiler(S, tw)            # theiler function in crqa package
    
    # assign 0 value to upper or lower part of RP
    if (side == "upper") {
      S = as.matrix(S)
      S[lower.tri(S, diag = TRUE)] = 0
      S = Matrix(S, sparse = TRUE)
    }
    if (side == "lower") {
      S = as.matrix(S)
      S[upper.tri(S, diag = TRUE)] = 0
      S = Matrix(S, sparse = TRUE)
    }
    if (side == "both") {
      S = S
    }
    
    # MY PIECE OF CODE FOR CATEGORICAL ENTROPY 
    # NB. this work appropriately only if S is a complete Recurrence Plot
    # (i.e. side == "both" above). Needs correction for the other situations.
    
    cat <- catEntropy(S)
    catentropy <- cat$Entropy
    bka <- cat$Blocks
    bka <- bka[bka[,5] > 1 & bka[,6] > 1, 7]
    distBlocks <- as.data.frame(table(bka))
    
    # extract non-zero diagonal lines in the form of a matrix B
    spdiagonalize = spdiags(S)   
    B = spdiagonalize$B
    
    # Computation of number of recurrence points and recurrence rate
    numrecurs = length(which(B == TRUE))
    percentrecurs = (numrecurs / ((v1l * v2l))) * 100
    
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
      numdiaglines = length(diaglines)
      maxline = max(diaglines)
      meanline = mean(diaglines)
      tabled = as.data.frame(table(diaglines))
      total = sum(tabled$Freq)
      p = tabled$Freq / total
      del = which(p == 0)
      if (length(del) > 0) {
        p = p[-del]
      }
      entropy = -sum(p * log(p))
      relEntropy = entropy / (-1 * log(1 / nrow(tabled)))
      pdeter = sum(diaglines) / numrecurs * 100
      restt = tt(S, minvertline, whiteline) # specific function in crqa
      lam = restt$lam
      TT = restt$TT
    }
    else {
      numdiaglines = 0
      maxline = 0
      pdeter = NA
      entropy = NA
      relEntropy = NA
      catENTR = catentropy
      distDLine = 0
      distBlocks = 0
      meanline = 0
      lam = 0
      TT = 0
      RP = NA
    }
    results = list(
      RR = percentrecurs,
      DET = pdeter,
      NRLINE = numdiaglines,
      maxL = maxline,
      L = meanline,
      ENTR = entropy,
      rENTR = relEntropy,
      catENTR = catentropy,
      distDLines = tabled,
      distBlocks = distBlocks,
      LAM = lam,
      TT = TT,
      RP = S
    )
  }
  
  # In the case of not a single recurrence point
  else {
    results = list(
      RR = 0,
      DET = NA,
      NRLINE = 0,
      maxL = 0,
      L = 0,
      ENTR = NA,
      rENTR = NA,
      catENTR = catentropy,
      distDLine = 0,
      distBlocks = 0,
      LAM = NA,
      TT = NA,
      RP = NA
    )
  }
  
  # What it spits out
  return(results)
}
