# test-minibatch.R
# author: Cristian Castiglione
# creation: 07/10/2023
# last change: 07/10/2023

## Workspace setup ----
rm(list = ls())
graphics.off()

# Package compilation and import
devtools::load_all()


r_get_next = function (iter, n, rnd) {
  idx = -1
  tovisit = seq(from = 0, to = n-1, by = 1)
  visited = c()
  if (iter > 0) {
    for (i in 1:iter) {
      if (length(tovisit) == 0) {
        tovisit = visited
        visited = c()
      }
      if (rnd) {
        idx = sample(tovisit, 1, replace = FALSE)
      } else {
        idx = tovisit[1]
      }
      j = which(tovisit == idx)
      tovisit = tovisit[-j]
      visited = c(visited, idx)
    }
  }
  list(idx = idx, tovisit = tovisit, visited = visited)
}


## Test: get_chunk() ----
sgdGMF::c_get_chunk(3, 10, 3, FALSE)

## Test: get_chunks() ----
sgdGMF::c_get_chunks(0:5, 10, 3, TRUE)

## Test: get_next() ----
print.pile = function (pile) {
  cat("idx =", pile$idx, "\n")
  cat("tovisit =", drop(pile$tovisit), "\n")
  cat("visited =", drop(pile$visited), "\n")
}

k = 0
print.pile(        r_get_next(k, 5, FALSE))
print.pile(sgdGMF::c_get_next(k, 5, FALSE))
k = k+1

## End of file ----
