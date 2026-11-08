# Parameters

base_params <- list(
  S0 = 100,
  K = 100,
  r = 0.05,
  sigma = 0.2,
  TE = 1
)

check_model <- function(S0, K, r, sigma, TE){
  pass <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  stopifnot(
    "S0 must be a single number > 0" = pass(S0) && S0 > 0,
    "K must be a single number > 0" = pass(K) && K > 0,
    "r must be a single finite number" = pass(r),
    "sigma must be a single number > 0" = pass(sigma) && sigma > 0,
    "TE must be a single number > 0" = pass(TE) && TE > 0
  )
}

################################################################################

# Point Generation

qmc_points <- function(
    type = c("halton", "ghalton", "sobol"),
    n,
    d,
    ghalton_seed = 67L, # fixed shift, NULL = new shift per call
    ...
){
  type <- match.arg(type)
  stopifnot(
    length(n) == 1, n >= 1, n == round(n),
    length(d) == 1, d >= 1, d == round(d)
  )
  
  U <- switch(
    type,
    halton = ghalton_points(
      n,
      d,
      seed = ghalton_seed,
      method = "halton"
    ),
    ghalton = ghalton_points(
      n,
      d,
      seed = ghalton_seed,
      method = "generalized"
    ),
    sobol = sobol_points(n, d, ...)
  )
  
  if(is.null(dim(U))) 
    U <- matrix(U, nrow = n, ncol = d)
  U
}

ghalton_points <- function(
    n,
    d,
    seed,
    method = c("halton", "generalized")
){
  method <- match.arg(method)
  stopifnot(
    "qrng::ghalton allows at most d = 360 dimensions" = d <= 360
  )
  if(is.null(seed))
    return(qrng::ghalton(n, d, method = method))
  withr::with_seed(
    seed,
    qrng::ghalton(n, d, method = method)
  )
}

sobol_points <- function(
    n,
    d,
    randomize = c(
      "none",
      "digital.shift",
      "Owen",
      "Faure.Tezuka",
      "Owen.Faure.Tezuka"
    ),
    ...
){
  randomize <- match.arg(randomize)
  stopifnot(
    n <= 2^31 - 1,
    d <= 16510
  )
  qrng::sobol(n, d, randomize = randomize, ...)
}

################################################################################

# CP-Rotation

cp_rotation <- function(
    P,
    delta = runif(ncol(P))
){
  stopifnot(
    is.matrix(P),
    length(delta) == ncol(P),
    "delta must lie in [0,1)" = all(delta >= 0 & delta < 1)
  )
  for(j in seq_len(ncol(P)))
    P[, j] <- (P[, j] + delta[j]) %% 1
  P[P == 0] <- .Machine$double.eps / 2 # 2^-53; machine epsilon
  # =  1 - largest double below 1;
  # since (p + delta) %% 1 can become exactly 0 in floating point arithmetic. 
  # prevents qnorm(0) = -inf
  P
}


################################################################################

# Paths

paths <- function(
    Z,
    t,
    construction = c("sequential", "bb")
){
  construction <- match.arg(construction)
  d <- length(t)
  stopifnot(
    is.matrix(Z),
    ncol(Z) == d,
    t[1] > 0,
    all(diff(t) > 0),
    "Z must be finite (qnorm(0) = -Inf)" = all(is.finite(Z))
  )
  tt <- c(0, t) # tt[i + 1] = t_i with t_0 = 0
  W <- matrix(0, nrow(Z), ncol = d)
  
  if(construction == "sequential"){
    W[, 1] <- sqrt(tt[2]) * Z[, 1]
    for(j in seq_len(d)[-1])
      W[, j] <- W[, j - 1] + sqrt(tt[j + 1] - tt[j]) * Z[, j]
  } else{
    W[, d] <- sqrt(tt[d + 1]) * Z[, 1]
    left <- 0L
    right <- d
    i <- 1L
    k <- 2L
    while(i <= length(left)){
      l <- left[i]
      r <- right[i]
      i <- i + 1L
      if(r - l < 2L) next
      m <- (l + r) %/% 2L
      
      tl <- tt[l + 1]
      tm <- tt[m + 1]
      tr <- tt[r + 1]
      
      Wl <- if(l == 0L)
        0 else
          W[, l]
      W[, m] <- ((tr - tm) * Wl + (tm - tl) * W[, r]) / (tr - tl) + 
        sqrt(((tm - tl) * (tr - tm)) / (tr - tl)) * Z[, k]
      k <- k + 1L
      left <- c(left, l, m)
      right <- c(right, m, r)
    }
  }
  W
}

################################################################################

# Payoff Function

payoff <- function(
    W,
    t,
    S0 = base_params$S0,
    K = base_params$K,
    r = base_params$r,
    sigma = base_params$sigma
){
  d <- length(t)
  stopifnot(is.matrix(W), ncol(W) == d)
  check_model(S0, K, r, sigma, t[d])
  sum_S <- numeric(nrow(W))
  sum_log <- numeric(nrow(W))
  for(j in seq_len(d)){
    logS <- log(S0) + (r - 0.5 * sigma^2) * t[j] + sigma * W[, j]
    sum_S <- sum_S + exp(logS)
    sum_log <- sum_log + logS
  }
  disc <- exp(-r * t[d])
  cbind(
    arithmetic = disc * pmax(sum_S / d - K, 0),
    geometric = disc * pmax(exp(sum_log / d) - K, 0)
  )
}

geometric_price <- function(
    d,
    S0 = base_params$S0,
    K = base_params$K,
    r = base_params$r,
    sigma = base_params$sigma,
    TE = base_params$TE
){
  stopifnot(
    "d must be a positive integer" = 
      is.numeric(d) && length(d) >= 1 && all(d >= 1) && all(d == round(d))
  )
  check_model(S0, K, r, sigma, TE)
  
  mu <- log(S0) + (r - 0.5 * sigma^2) * ((d + 1) * TE / (2 * d)) # mu_d
  v2 <- (sigma^2 * TE * (d + 1) * (2 * d + 1)) / (6 * d^2) # v_d^2
  v <- sqrt(v2)
  d1 <- (mu - log(K) + v2) / v
  d2 <- d1 - v
  exp(-r * TE) * (exp(mu + 0.5 * v2) * pnorm(d1) - K * pnorm(d2))
}

################################################################################

# Unexplained share

# shares(s)
# = (sigma_f^2 - tau_s^2) / sigma_f^2 = E[f(U) - f(U^(s)))^2] / (2 sigma_f^2)

unexplained_share <- function(
    d,
    construction,
    s = 1:d,
    n = 2^16,
    which_payoff = "arithmetic",
    seed = 1,
    params = base_params
){
  stopifnot(
    "d must be a single positive integer" =
      is.numeric(d) && length(d) == 1 && d == round(d),
    "s must be  coordinates in 1,..., d" =
      all(s %in% seq_len(d))
  )
  
  set.seed(seed)
  t <- (1:d) * params$TE / d
  f <- function(U)
    payoff(
      paths(qnorm(U), t, construction),
      t,
      params$S0,
      params$K,
      params$r,
      params$sigma
    )[, which_payoff]
  
  U <- matrix(runif(n * d), n, d)
  V <- matrix(runif(n * d), n, d)
  
  f_U <- f(U)
  sigma2_f <- var(f_U)
  
  share <- numeric(length(s))
  for(i in seq_along(s)){
    k <- s[i]
    U_s <- V
    U_s[, 1:k] <- U[, 1:k]
    share[i] <- 0.5 * mean((f_U - f(U_s))^2) / sigma2_f
  }
  data.frame(
    d = d,
    construction = factor(
      construction,
      levels = c("sequential", "bb"),
      labels = c("Sequential", "Brownian Bridge")
    ),
    s = s,
    share = share
  )
}








