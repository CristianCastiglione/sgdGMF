
# Sample the minibatch partition for all the iterations of the algorithm
sample.minibatch = function (nobs, size, randomize = TRUE) {
  # do you want to shuffle the data or not?
  if (randomize) {
    # if the data are ordered, it might be useful to shuffle the observations
    idx = sample.int(nobs, nobs, replace = FALSE)
  } else {
    # if the data are already randomized, we can keep them in the original order
    idx = 1:nobs
  }
  # divide the data in chunks of (almost) the same size
  chunks = split(idx, ceiling(seq_along(idx) / size))
  # return the obtained chunks
  return (chunks)
}

# Select the minibatch to use at each iteration
select.minibatch = function (iter, nmb) {
  mod = iter %% nmb
  idx = if (mod == 0) nmb else mod
  return (idx)
}

# Sample the minibatch indeces at each iteration
## sample.minibatch.old = function (n, size) {
##   # ngroup = floor(n / size)
##   # labels = cut(1:n, ngroup, labels = FALSE)
##   # groups = split(1:n, idx)
##   # return (groups)
##   sample(1:n, size = size, replace = FALSE)
## }

# Select the minibatch to use among the ones collected in set
# Return the choosen minibatch as idx, the remaining set as
# setnew = setnew \ {i} and the complementary set as setold = {setold, i}
no.replace.selection = function (setnew, setold) {
  if (length(setnew) == 0) {setnew = setold}
  idx = sample(setnew, 1, replace = FALSE)
  setnew = setdiff(setnew, idx)
  setold = c(setold, idx)
  list(idx = idx, setnew = setnew, setold = setold)
}

# Sample the minibatch indeces based on a group variable
stratified.minibatch = function (n, size, group) {
  # ngroup = floor(n / size)
  # labels = cut(1:n, ngroup, labels = FALSE)
  # groups = split(1:n, idx)
  # return (groups)
  sample(1:n, size = size, replace = FALSE)
}

# Compute the linear predictor
## lin.pred.gmf = function (U, V,
##                          X = NULL, Z = NULL, Bx = NULL, Bz = NULL,
##                          idxr = NULL, idxc = NULL, trasp = FALSE) {
##
##   eta = NULL
##   if (is.null(idxr) & is.null(idxc)) {
##     eta = tcrossprod(U, V)
##     if (!is.null(X)) eta = eta + tcrossprod(X, Bx)
##     if (!is.null(Z)) eta = eta + tcrossprod(Bz, Z)
##   }
##   if (!is.null(idxr) & is.null(idxc)) {
##     eta = tcrossprod(U[idxr, , drop = FALSE], V)
##     if (!is.null(X)) eta = eta + tcrossprod(X[idxr, , drop = FALSE], Bx)
##     if (!is.null(Z)) eta = eta + tcrossprod(Bz[idxr, , drop = FALSE], Z)
##   }
##   if (is.null(idxr) & !is.null(idxc)) {
##     eta = tcrossprod(U, V[idxc, , drop = FALSE])
##     if (!is.null(X)) eta = eta + tcrossprod(X, Bx[idxc, , drop = FALSE])
##     if (!is.null(Z)) eta = eta + tcrossprod(Bz, Z[idxc, , drop = FALSE])
##   }
##   if (!is.null(idxr) & !is.null(idxc)) {
##     eta = tcrossprod(U[idxr, , drop = FALSE], V[idxc, , drop = FALSE])
##     if (!is.null(X)) eta = eta + tcrossprod(X[idxr, , drop = FALSE], Bx[idxc, , drop = FALSE])
##     if (!is.null(Z)) eta = eta + tcrossprod(Bz[idxr, , drop = FALSE], Z[idxc, , drop = FALSE])
##   }
##
##   if (trasp) eta = t(eta)
##
##   return (eta)
## }

# Compute the GLM residuals
residuals.gmf = function (Y, eta, family, type = "pearson") {

  if (!(type %in% c("raw", "pearson", "deviance"))) {
    stop("Residual type not allowed.")
  }

  n = nrow(Y)
  m = ncol(Y)
  res = matrix(NA, nrow = n, ncol = m)
  mu = matrix(family$linkinv(eta), nrow = n, ncol = m)

  if (type == "raw") {
    res = Y - mu
  }
  if (type == "pearson") {
    var = matrix(family$variance(mu), nrow = n, ncol = m)
    res = (Y - mu) / var
  }
  if (type == "deviance") {
    dev = pointwise.deviance(mu, Y, family)
    res = sign(Y - mu) * sqrt(abs(dev))
  }

  return (res)
}

# Update the rate parameter following the Robbins-Monro rule
update.sgd.rate = function (iter, rate, decay) {
  return (rate / (1 + decay * rate * iter)^.75)
}

# Calculate the first and second otder derivatives of the log-likelihood
# function with respect to the linear predictor (eta)
get.eta.deriv = function (Y, eta, family, trasp = FALSE) {

  n = nrow(Y)
  m = ncol(Y)

  # Calcualte the mean, variance and inverse link derivative
  mu = matrix(family$linkinv(eta), nrow = n, ncol = m)
  var = matrix(family$variance(mu), nrow = n, ncol = m)
  dginv = matrix(family$mu.eta(eta), nrow = n, ncol = m)
  ratio = matrix(dginv / var, nrow = n, ncol = m)

  # Calculate the first and second order derivatives of the
  # unpenalized log-likelihood wrt the linear predictor
  deta = ratio * (Y - mu)
  ddeta = ratio * dginv

  if (trasp) {
    deta = t(deta)
    ddeta = t(ddeta)
  }

  # Return the derivatives
  list(deta = deta, ddeta = ddeta)
}

get.dispersion = function (y, mu, family) {
  colMeans((y - mu)^2 / family$variance(mu))
}

# Compute the smoothed gradients for the patameters
update.sgd.grad = function (U, V, deta, gradt, rate = 0.9,
                            scale = 1, penalty = 0,
                            trasp = FALSE) {

  # compute the cuttent gradient
  grad = - deta %*% V + penalty * U

  # smooth the historical and current gradients
  grad = (1 - rate) * gradt + rate * grad

  # transpose the result if required
  if (trasp) grad = t(grad)

  # return the smoothed gradient
  return (grad)
}

# Compute the smoothed hessians for the parameters
update.sgd.hess = function (U, V, ddeta, hesst, rate = 0.9,
                            scale = 1, penalty = 0, damping = 0,
                            trasp = FALSE) {

  # compute the cuttent Hessian
  hess = ddeta %*% V**2 + penalty + damping

  # smooth the historical and current Hessians
  hess = (1 - rate) * hesst + rate * hess

  # transpose the result if required
  if (trasp) hess = t(hess)

  # return the smoothed gradient
  return (hess)
}

# Compute the parameter updates
update.sgd.params = function (Ut, gradt, hesst, ratet = 0.1) {
  Ut - ratet * (gradt / hesst)
}

# Smooth the parameter estimates
smooth.sgd.params = function (U, Ut, burn, iter) {

  # If the warmup period is concluded, average the estimated parametrs
  # for all the subsequent iterations using a uniform weigthing
  if (iter > burn) {
    ratet = 1 / (iter - burn)
    U = (1 - ratet) * Ut + ratet * U
  }

  return (U)
}

# Print the optimization status
print.status = function (iter, deviance, change, exetime,
                         rowfrac = NULL, colfrac = NULL) {

  if (iter == 1) {
    if (is.null(rowfrac) & is.null(colfrac)) {
      cat(" Iteration   ",
          " Deviance  ",
          " Change  ",
          " Exe.Time  ",
          "\n", sep = "")
    } else {
      cat(" Iteration   ",
          " Deviance  ",
          " Change  ",
          " Exe.Time  ",
          " Scanned.rows ",
          " Scanned.cols ",
          "\n", sep = "")
    }
  }

  if (exetime / 60 < 1) {
    if (is.null(rowfrac) | is.null(colfrac)) {
      cat(paste(gettextf(" %9d ", iter),
                gettextf(" %9.4f ", deviance),
                gettextf(" %2.4f ", change),
                gettextf(" %6.2f s", exetime),
                "\n", sep = " "))
    } else {
      cat(paste(gettextf(" %9d ", iter),
                gettextf(" %9.4f ", deviance),
                gettextf(" %2.4f ", change),
                gettextf(" %6.2f s", exetime),
                gettextf(" %9d/100", rowfrac),
                gettextf(" %8d/100", colfrac),
                "\n", sep = " "))
    }
  } else {
    if (is.null(rowfrac) | is.null(colfrac)) {
      cat(paste(gettextf(" %9d ", iter),
                gettextf(" %9.4f ", deviance),
                gettextf(" %2.4f ", change),
                gettextf(" %6.2f m", exetime / 60),
                gettextf(" %9d/100", rowfrac),
                gettextf(" %8d/100", colfrac),
                "\n", sep = " "))
    } else {
      cat(paste(gettextf(" %9d ", iter),
                gettextf(" %9.4f ", deviance),
                gettextf(" %2.4f ", change),
                gettextf(" %6.2f m", exetime / 60),
                gettextf(" %9d/100", rowfrac),
                gettextf(" %8d/100", colfrac),
                "\n", sep = " "))
    }
  }
}


# Print the optimization status
print.status.2 = function (iter, deviance, change, exetime, scanned = NULL) {

  if (iter == 1) {
    cat(" Iteration   ",
        " Deviance  ",
        " Change  ",
        " Exe.Time  ",
        " Scanned.data ",
        "\n", sep = "")
  }

  if (exetime / 60 < 1) {
    if (is.null(scanned)) {
      cat(paste(gettextf(" %9d ", iter),
                gettextf(" %9.4f ", deviance),
                gettextf(" %2.4f ", change),
                gettextf(" %6.2f s", exetime),
                "\n", sep = " "))
    } else {
      cat(paste(gettextf(" %9d ", iter),
                gettextf(" %9.4f ", deviance),
                gettextf(" %2.4f ", change),
                gettextf(" %6.2f s", exetime),
                gettextf(" %9d/100", scanned),
                "\n", sep = " "))
    }
  } else {
    if (is.null(scanned)) {
      cat(paste(gettextf(" %9d ", iter),
                gettextf(" %9.4f ", deviance),
                gettextf(" %2.4f ", change),
                gettextf(" %6.2f m", exetime / 60),
                "\n", sep = " "))
    } else {
      cat(paste(gettextf(" %9d ", iter),
                gettextf(" %9.4f ", deviance),
                gettextf(" %2.4f ", change),
                gettextf(" %6.2f m", exetime / 60),
                gettextf(" %9d/100", scanned),
                "\n", sep = " "))
    }
  }
}

progress.bar = function (progress, time, length = 50) {
  if (time > 60) {
    time = paste(floor(time / 60), "m", sep = "")
  } else {
    time = paste(floor(time), "s", sep = "")
  }
  if (progress < length) {
    bar = c(rep("=", progress), ">", rep(".", length - progress - 1))
    cat("\r |", bar, "| ", 2 * progress, "%, ", time, "  ", sep = "")
  } else {
    cat("\r |", rep("=", progress), "| 100%, ", time, "\n", sep = "")
  }
}
