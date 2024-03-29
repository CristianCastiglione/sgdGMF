# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

cpp.airwls.glmstep <- function(beta, y, X, familyname, linkname, offset, penalty) {
    .Call(`_sgdGMF_cpp_airwls_glmstep`, beta, y, X, familyname, linkname, offset, penalty)
}

cpp.airwls.glmfit <- function(beta, y, X, familyname, linkname, offset, penalty, nsteps = 100L, stepsize = 0.1, print = FALSE) {
    .Call(`_sgdGMF_cpp_airwls_glmfit`, beta, y, X, familyname, linkname, offset, penalty, nsteps, stepsize, print)
}

cpp.airwls.update <- function(beta, Y, X, familyname, linkname, idx, offset, penalty, transp = FALSE, nsteps = 100L, stepsize = 0.1, print = FALSE, parallel = FALSE, nthreads = 1L) {
    .Call(`_sgdGMF_cpp_airwls_update`, beta, Y, X, familyname, linkname, idx, offset, penalty, transp, nsteps, stepsize, print, parallel, nthreads)
}

cpp.fit.airwls <- function(Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter = 500L, nsteps = 1L, stepsize = 0.1, eps = 1e-08, nafill = 1L, tol = 1e-05, damping = 1e-03, verbose = TRUE, frequency = 10L, parallel = FALSE, nthreads = 1L) {
    .Call(`_sgdGMF_cpp_fit_airwls`, Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter, nsteps, stepsize, eps, nafill, tol, damping, verbose, frequency, parallel, nthreads)
}

cpp.fit.newton <- function(Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter = 500L, stepsize = 0.1, eps = 1e-08, nafill = 1L, tol = 1e-05, damping = 1e-03, verbose = TRUE, frequency = 10L, parallel = FALSE, nthreads = 1L) {
    .Call(`_sgdGMF_cpp_fit_newton`, Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter, stepsize, eps, nafill, tol, damping, verbose, frequency, parallel, nthreads)
}

cpp.fit.msgd <- function(Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter = 1000L, eps = 0.01, nafill = 10L, tol = 1e-08, size = 100L, burn = 0.75, rate0 = 0.01, decay = 0.01, damping = 1e-03, rate1 = 0.95, rate2 = 0.99, parallel = FALSE, nthreads = 1L, verbose = TRUE, frequency = 250L, progress = FALSE) {
    .Call(`_sgdGMF_cpp_fit_msgd`, Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter, eps, nafill, tol, size, burn, rate0, decay, damping, rate1, rate2, parallel, nthreads, verbose, frequency, progress)
}

cpp.fit.csgd <- function(Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter = 1000L, eps = 0.01, nafill = 10L, tol = 1e-08, size1 = 100L, size2 = 100L, burn = 0.75, rate0 = 0.01, decay = 0.01, damping = 1e-03, rate1 = 0.95, rate2 = 0.99, parallel = FALSE, nthreads = 1L, verbose = TRUE, frequency = 250L, progress = FALSE) {
    .Call(`_sgdGMF_cpp_fit_csgd`, Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter, eps, nafill, tol, size1, size2, burn, rate0, decay, damping, rate1, rate2, parallel, nthreads, verbose, frequency, progress)
}

cpp.fit.rsgd <- function(Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter = 1000L, eps = 0.01, nafill = 10L, tol = 1e-08, size1 = 100L, size2 = 100L, burn = 0.75, rate0 = 0.01, decay = 0.01, damping = 1e-03, rate1 = 0.95, rate2 = 0.99, parallel = FALSE, nthreads = 1L, verbose = TRUE, frequency = 250L, progress = FALSE) {
    .Call(`_sgdGMF_cpp_fit_rsgd`, Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter, eps, nafill, tol, size1, size2, burn, rate0, decay, damping, rate1, rate2, parallel, nthreads, verbose, frequency, progress)
}

cpp.fit.bsgd <- function(Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter = 1000L, eps = 0.01, nafill = 10L, tol = 1e-08, size1 = 100L, size2 = 100L, burn = 0.75, rate0 = 0.01, decay = 0.01, damping = 1e-03, rate1 = 0.95, rate2 = 0.99, parallel = FALSE, nthreads = 1L, verbose = TRUE, frequency = 250L, progress = FALSE) {
    .Call(`_sgdGMF_cpp_fit_bsgd`, Y, X, B, A, Z, U, V, familyname, linkname, ncomp, lambda, maxiter, eps, nafill, tol, size1, size2, burn, rate0, decay, damping, rate1, rate2, parallel, nthreads, verbose, frequency, progress)
}

cpp.deviance <- function(y, mu, familyname) {
    .Call(`_sgdGMF_cpp_deviance`, y, mu, familyname)
}

cpp.penalty <- function(u, p) {
    .Call(`_sgdGMF_cpp_penalty`, u, p)
}

pcc.gaussian.variance <- function(mu) {
    .Call(`_sgdGMF_cpp_gaussian_variance`, mu)
}

pcc.gaussian.initialize <- function(y) {
    .Call(`_sgdGMF_cpp_gaussian_initialize`, y)
}

pcc.gaussian.devresid <- function(y, mu) {
    .Call(`_sgdGMF_cpp_gaussian_devresid`, y, mu)
}

pcc.binomial.variance <- function(mu) {
    .Call(`_sgdGMF_cpp_binomial_variance`, mu)
}

pcc.binomial.initialize <- function(y) {
    .Call(`_sgdGMF_cpp_binomial_initialize`, y)
}

pcc.binomial.devresid <- function(y, mu) {
    .Call(`_sgdGMF_cpp_binomial_devresid`, y, mu)
}

pcc.poisson.variance <- function(mu) {
    .Call(`_sgdGMF_cpp_poisson_variance`, mu)
}

pcc.poisson.initialize <- function(y) {
    .Call(`_sgdGMF_cpp_poisson_initialize`, y)
}

pcc.poisson.devresid <- function(y, mu) {
    .Call(`_sgdGMF_cpp_poisson_devresid`, y, mu)
}

pcc.gamma.variance <- function(mu) {
    .Call(`_sgdGMF_cpp_gamma_variance`, mu)
}

pcc.gamma.initialize <- function(y) {
    .Call(`_sgdGMF_cpp_gamma_initialize`, y)
}

pcc.gamma.devresid <- function(y, mu) {
    .Call(`_sgdGMF_cpp_gamma_devresid`, y, mu)
}

pcc.negbinom.variance <- function(mu) {
    .Call(`_sgdGMF_cpp_negbinom_variance`, mu)
}

pcc.negbinom.initialize <- function(y) {
    .Call(`_sgdGMF_cpp_negbinom_initialize`, y)
}

pcc.negbinom.devresid <- function(y, mu) {
    .Call(`_sgdGMF_cpp_negbinom_devresid`, y, mu)
}

cpp.link.identity.linkfun <- function(mu) {
    .Call(`_sgdGMF_cpp_link_identity_linkfun`, mu)
}

cpp.link.identity.linkinv <- function(eta) {
    .Call(`_sgdGMF_cpp_link_identity_linkinv`, eta)
}

cpp.link.identity.mueta <- function(eta) {
    .Call(`_sgdGMF_cpp_link_identity_mueta`, eta)
}

cpp.link.logit.linkfun <- function(mu) {
    .Call(`_sgdGMF_cpp_link_logit_linkfun`, mu)
}

cpp.link.logit.linkinv <- function(eta) {
    .Call(`_sgdGMF_cpp_link_logit_linkinv`, eta)
}

cpp.link.logit.mueta <- function(eta) {
    .Call(`_sgdGMF_cpp_link_logit_mueta`, eta)
}

cpp.link.probit.linkfun <- function(mu) {
    .Call(`_sgdGMF_cpp_link_probit_linkfun`, mu)
}

cpp.link.probit.linkinv <- function(eta) {
    .Call(`_sgdGMF_cpp_link_probit_linkinv`, eta)
}

cpp.link.probit.mueta <- function(eta) {
    .Call(`_sgdGMF_cpp_link_probit_mueta`, eta)
}

cpp.link.cauchy.linkfun <- function(mu) {
    .Call(`_sgdGMF_cpp_link_cauchy_linkfun`, mu)
}

cpp.link.cauchy.linkinv <- function(eta) {
    .Call(`_sgdGMF_cpp_link_cauchy_linkinv`, eta)
}

cpp.link.cauchy.mueta <- function(eta) {
    .Call(`_sgdGMF_cpp_link_cauchy_mueta`, eta)
}

cpp.link.cloglog.linkfun <- function(mu) {
    .Call(`_sgdGMF_cpp_link_cloglog_linkfun`, mu)
}

cpp.link.cloglog.linkinv <- function(eta) {
    .Call(`_sgdGMF_cpp_link_cloglog_linkinv`, eta)
}

cpp.link.cloglog.mueta <- function(eta) {
    .Call(`_sgdGMF_cpp_link_cloglog_mueta`, eta)
}

cpp.link.log.linkfun <- function(mu) {
    .Call(`_sgdGMF_cpp_link_log_linkfun`, mu)
}

cpp.link.log.linkinv <- function(eta) {
    .Call(`_sgdGMF_cpp_link_log_linkinv`, eta)
}

cpp.link.log.mueta <- function(eta) {
    .Call(`_sgdGMF_cpp_link_log_mueta`, eta)
}

cpp.link.inverse.linkfun <- function(mu) {
    .Call(`_sgdGMF_cpp_link_inverse_linkfun`, mu)
}

cpp.link.inverse.linkinv <- function(eta) {
    .Call(`_sgdGMF_cpp_link_inverse_linkinv`, eta)
}

cpp.link.inverse.mueta <- function(eta) {
    .Call(`_sgdGMF_cpp_link_inverse_mueta`, eta)
}

cpp.link.sqrt.linkfun <- function(mu) {
    .Call(`_sgdGMF_cpp_link_sqrt_linkfun`, mu)
}

cpp.link.sqrt.linkinv <- function(eta) {
    .Call(`_sgdGMF_cpp_link_sqrt_linkinv`, eta)
}

cpp.link.sqrt.mueta <- function(eta) {
    .Call(`_sgdGMF_cpp_link_sqrt_mueta`, eta)
}

cpp.get.chunk <- function(iter, n, size, randomize) {
    .Call(`_sgdGMF_cpp_get_chunk`, iter, n, size, randomize)
}

cpp.get.chunks <- function(iters, n, size, randomize) {
    .Call(`_sgdGMF_cpp_get_chunks`, iters, n, size, randomize)
}

cpp.get.next <- function(iter, n, rnd) {
    .Call(`_sgdGMF_cpp_get_next`, iter, n, rnd)
}

cpp.make.link.family <- function(familyname, linkname) {
    invisible(.Call(`_sgdGMF_cpp_make_link_family`, familyname, linkname))
}

cpp.get.data.bounds <- function(eps, ymin, ymax, familyname, linkname) {
    .Call(`_sgdGMF_cpp_get_data_bounds`, eps, ymin, ymax, familyname, linkname)
}

cpp.get.uv.penalty <- function(pen, p, q, d) {
    .Call(`_sgdGMF_cpp_get_uv_penalty`, pen, p, q, d)
}

cpp.get.uv.indices <- function(p, q, d) {
    .Call(`_sgdGMF_cpp_get_uv_indices`, p, q, d)
}

cpp.sample.minibatch <- function(n, size, randomize) {
    .Call(`_sgdGMF_cpp_sample_minibatch`, n, size, randomize)
}

cpp.select.minibatch <- function(iter, nchunks) {
    .Call(`_sgdGMF_cpp_select_minibatch`, iter, nchunks)
}

cpp.dabsmax <- function(u, v) {
    .Call(`_sgdGMF_cpp_dabsmax`, u, v)
}

cpp.vabsmax <- function(u, v) {
    .Call(`_sgdGMF_cpp_vabsmax`, u, v)
}

cpp.trim <- function(x, a, b) {
    .Call(`_sgdGMF_cpp_trim`, x, a, b)
}

cpp.xlogx <- function(x) {
    .Call(`_sgdGMF_cpp_xlogx`, x)
}

cpp.log1pexp <- function(x) {
    .Call(`_sgdGMF_cpp_log1pexp`, x)
}

cpp.log1mexp <- function(x) {
    .Call(`_sgdGMF_cpp_log1mexp`, x)
}

cpp.logit <- function(x) {
    .Call(`_sgdGMF_cpp_logit`, x)
}

cpp.expit <- function(x) {
    .Call(`_sgdGMF_cpp_expit`, x)
}

cpp.expit2 <- function(x) {
    .Call(`_sgdGMF_cpp_expit2`, x)
}

cpp.expitn <- function(x, n = 1) {
    .Call(`_sgdGMF_cpp_expitn`, x, n)
}

cpp.cloglog <- function(x) {
    .Call(`_sgdGMF_cpp_cloglog`, x)
}

cpp.cexpexp <- function(x) {
    .Call(`_sgdGMF_cpp_cexpexp`, x)
}

cpp.loglog <- function(x) {
    .Call(`_sgdGMF_cpp_loglog`, x)
}

cpp.expexp <- function(x) {
    .Call(`_sgdGMF_cpp_expexp`, x)
}

cpp.pdfn <- function(x) {
    .Call(`_sgdGMF_cpp_pdfn`, x)
}

cpp.cdfn <- function(x) {
    .Call(`_sgdGMF_cpp_cdfn`, x)
}

cpp.logpdfn <- function(x) {
    .Call(`_sgdGMF_cpp_logpdfn`, x)
}

cpp.logcdfn <- function(x) {
    .Call(`_sgdGMF_cpp_logcdfn`, x)
}

cpp.gamma <- function(x) {
    .Call(`_sgdGMF_cpp_gamma`, x)
}

cpp.loggamma <- function(x) {
    .Call(`_sgdGMF_cpp_loggamma`, x)
}

cpp.digamma <- function(x) {
    .Call(`_sgdGMF_cpp_digamma`, x)
}

cpp.trigamma <- function(x) {
    .Call(`_sgdGMF_cpp_trigamma`, x)
}

cpp.beta <- function(x, y) {
    .Call(`_sgdGMF_cpp_beta`, x, y)
}

cpp.logbeta <- function(x, y) {
    .Call(`_sgdGMF_cpp_logbeta`, x, y)
}

cpp.dibeta <- function(x, y) {
    .Call(`_sgdGMF_cpp_dibeta`, x, y)
}

cpp.tribeta <- function(x, y) {
    .Call(`_sgdGMF_cpp_tribeta`, x, y)
}

cpp.hinge <- function(x) {
    .Call(`_sgdGMF_cpp_hinge`, x)
}

cpp.dirac <- function(x, a = 0) {
    .Call(`_sgdGMF_cpp_dirac`, x, a)
}

cpp.step <- function(x, a = 0, lower = TRUE) {
    .Call(`_sgdGMF_cpp_step`, x, a, lower)
}

cpp.vech <- function(A) {
    .Call(`_sgdGMF_cpp_vech`, A)
}

