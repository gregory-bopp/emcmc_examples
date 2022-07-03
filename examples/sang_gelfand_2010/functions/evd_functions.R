qgev <- function(p,
                  loc = 0,
                  scale = 1,
                  shape = 0) {
  if (length(loc) > 1) {
    loc <- as.matrix(loc)
    if (ncol(p) != ncol(loc)) {
      loc <- matrix(loc, nrow = nrow(p), ncol = ncol(p))
    }
    out <- matrix(NA, nrow = nrow(p), ncol = ncol(p))
    for (i in 1:nrow(out)) {
      for (j in 1:ncol(out)) {
        out[i, j] <-
          qevd(p[i, j],
               loc = loc[i, j],
               scale = scale,
               shape = shape)
      }
    }
    return(out)
  }
  else{
    return(qevd(
      p,
      loc = loc,
      scale = scale,
      shape = shape
    ))
  }
}

pgev <- function(q,
                  loc = 0,
                  scale = 1,
                  shape = 0) {
  if (length(loc) > 1) {
    loc <- as.matrix(loc)
    if (ncol(q) != ncol(loc)) {
      loc <- matrix(loc, nrow = nrow(q), ncol = ncol(q))
    }
    out <- matrix(NA, nrow = nrow(q), ncol = ncol(q))
    for (i in 1:nrow(out)) {
      for (j in 1:ncol(out)) {
        out[i, j] <-
          pevd(q[i, j],
               loc = loc[i, j],
               scale = scale,
               shape = shape)
      }
    }
    return(out)
  }
  else{
    return(pevd(
      q,
      loc = loc,
      scale = scale,
      shape = shape
    ))
  }
}


dgev <- function(x,
                 loc = 0,
                 scale = 1,
                 shape = 0,
                 log = FALSE) {
  if (length(loc) > 1) {
    loc <- as.matrix(loc)
    if (ncol(x) != ncol(loc)) {
      loc <- matrix(loc, nrow = nrow(x), ncol = ncol(x))
    }
    out <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (i in 1:nrow(out)) {
      for (j in 1:ncol(out)) {
        out[i, j] <-
          devd(x[i, j],
               loc = loc[i, j],
               scale = scale,
               shape = shape,
               log = log)
      }
    }
    return(out)
  }
  else{
    return(pevd(
      x,
      loc = loc,
      scale = scale,
      shape = shape
    ))
  }
}
