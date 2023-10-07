# test-misc.R
# author: Cristian Castiglione
# creation: 02/10/2023
# last change: 06/10/2023

## Workspace setup ----
rm(list = ls())
graphics.off()

# Package compilation and import
devtools::load_all()


## Test: get_data_bounds() ----
{
  ymin = 0; ymax = 1; eps = 0.01
  r.bounds = binomial(link = "probit")$linkfun(c(ymin+eps*(ymax-ymin), ymax-eps*(ymax-ymin)))
  c.bounds = drop(sgdGMF::c_get_data_bounds(eps, ymin, ymax, "binomial", "probit")$etalim)
  print(all.equal(r.bounds, c.bounds))
}

{
  ymin = 0; ymax = 1; eps = 0.01
  r.bounds = binomial(link = "logit")$linkfun(c(ymin+eps*(ymax-ymin), ymax-eps*(ymax-ymin)))
  c.bounds = drop(sgdGMF::c_get_data_bounds(eps, ymin, ymax, "binomial", "logit")$etalim)
  print(all.equal(r.bounds, c.bounds))
}

## Test: get_uv_indices() ----
{
  p = 3; q = 1; d = 2
  r.idx = list(idu = c(p:(p+q-1), (p+q):(p+q+d-1)), idv = c(0:(p-1), (p+q):(p+q+d-1)))
  c.idx = sgdGMF::c_get_uv_indices(p, q, d)
  print(all.equal(r.idx$idu, drop(c.idx$idu)))
  print(all.equal(r.idx$idv, drop(c.idx$idv)))
}

## Test: get_uv_penalty() ----
{
  p = 3; q = 1; d = 2; pen = c(1:4)
  r.pen = list(penu = c(rep(0,p), rep(pen[1],q), rep(pen[3],d)),
               penv = c(rep(pen[2],p), rep(0,q), rep(pen[4],d)))
  c.pen = sgdGMF::c_get_uv_penalty(pen, p, q, d)
  print(all.equal(r.pen$penu, drop(c.pen$penu)))
  print(all.equal(r.pen$penv, drop(c.pen$penv)))
}

## Test: sample_minibatch() ----
{
  n = 9; size = 3; randomize = FALSE
  r.chunks = sgdGMF::sample.minibatch(n, size, randomize)
  c.chunks = sgdGMF::c_sample_minibatch(n, size, randomize)

  flag = TRUE
  for (h in 1:ceiling(n / size)) {
    flagh = all.equal(r.chunks[[h]], c.chunks[[h]]+1)
    flag = flag && flagh
  }
  print(flag)
}


{
  n = 11; size = 3; randomize = FALSE
  r.chunks = sgdGMF::sample.minibatch(n, size, randomize)
  c.chunks = sgdGMF::c_sample_minibatch(n, size, randomize)

  flag = TRUE
  for (h in 1:ceiling(n / size)) {
    flagh = all.equal(r.chunks[[h]], c.chunks[[h]]+1)
    flag = flag && flagh
  }
  print(flag)
}

## Test: select_chunk ----
{
  iter = 10; nchunks = 3
  r.idx = sgdGMF::select.minibatch(iter, nchunks)
  c.idx = sgdGMF::c_select_minibatch(iter, nchunks)
  print(all.equal(r.idx, c.idx+1))
}

## Test: get_chunks ----
{
  n = 9; size = 3; randomize = FALSE
  r.chunks = sgdGMF::sample.minibatch(n, size, randomize)
  c.chunks = sgdGMF::c_get_chunks(0:2, n, size, randomize)

  flag = TRUE
  for (h in 1:ceiling(n / size)) {
    flagh = all.equal(r.chunks[[h]], c.chunks[[h]]+1)
    flag = flag && flagh
  }
  print(flag)
}

{
  n = 10; size = 3; randomize = FALSE
  r.chunks = sgdGMF::sample.minibatch(n, size, randomize)
  c.chunks = sgdGMF::c_get_chunks(0:3, n, size, randomize)

  flag = TRUE
  for (h in 1:ceiling(n / size)) {
    flagh = all.equal(r.chunks[[h]], c.chunks[[h]]+1)
    flag = flag && flagh
  }
  print(flag)
}

{
  n = 11; size = 3; randomize = FALSE
  r.chunks = sgdGMF::sample.minibatch(n, size, randomize)
  c.chunks = sgdGMF::c_get_chunks(0:8, n, size, randomize)

  flag = TRUE
  for (h in 1:ceiling(n / size)) {
    flagh = all.equal(r.chunks[[h]], c.chunks[[h]]+1)
    flag = flag && flagh
  }
  print(flag)
}

## End of file ----

