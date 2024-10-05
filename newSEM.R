new_NullSEM <- function(x = double(), convergence = FALSE) {
  stopifnot(is.double(x))
  structure(x,
    class = "NullSEM",
    converged = convergence
  )
}

inspect.NullSEM <- function(x, what) {
  attr(x, what)
}
