source("R_scripts/functions.R")

################################################################################

# Setup

dims <- c(1, 4, 12, 52, 252)
ns <- 2^(8:16) # sample sizes
M <- 100 # random shifts

dir.create("results", showWarnings = FALSE)

################################################################################

run_experiment <- function(
    method,
    d,
    construction,
    ns,
    M,
    params = base_params,
    shift = TRUE
){
  method <- match.arg(method, c("mc", "halton", "ghalton", "sobol"))
  construction <- match.arg(construction, c("sequential", "bb"))
  stopifnot(
    "d must be a single positive integer" =
      is.numeric(d) && length(d) == 1 && d >= 1 && d == round(d),
    "ns must be positive integers" =
      is.numeric(ns) && length(ns) >= 1 && all(ns >= 1) && all(ns == round(ns)),
    "M must be an integer >= 2" =
      is.numeric(M) && length(M) == 1 && M >= 2 && M == round(M),
    "shift must be TRUE or FALSE" =
      is.logical(shift) && length(shift) == 1 && !is.na(shift)
    )
  check_model(
    params$S0,
    params$K,
    params$r,
    params$sigma,
    params$TE
  )
  
  t <- (1:d) * params$TE / d
  n_max <- max(ns)
  reps <- if(method != "mc" && !shift) 1 else M
  estimator <- if(method == "mc") "MC" else if(shift) "RQMC" else "QMC"

  set.seed(1000L * match(method, c("mc", "halton", "ghalton", "sobol")) + d)
  
  if(method != "mc")
    P <- qmc_points(method, n_max, d)
  
  out <- vector("list", reps)
  for(i in seq_len(reps)){
    if(method == "mc"){
      Z <- matrix(rnorm(n_max * d), n_max, d)
    } else if(shift){
      Z <- qnorm(cp_rotation(P))
    } else{
      Z <- qnorm(P)
    }
    
    W <- paths(Z, t, construction)
    X <- payoff(
      W,
      t,
      params$S0,
      params$K,
      params$r,
      params$sigma
    )
    
    out[[i]] <- data.frame(
      method = method,
      estimator = estimator,
      d = d,
      construction = construction,
      rep = i,
      n = ns,
      arithmetic = cumsum(X[, "arithmetic"])[ns] / ns,
      geometric = cumsum(X[, "geometric"])[ns] / ns
    )
  }
  do.call(rbind, out)
}

################################################################################

designs <- expand.grid(
  method = c("mc", "halton", "ghalton", "sobol"),
  d = dims,
  construction = c("sequential", "bb"),
  stringsAsFactors = FALSE
)

results <- list()
t0 <- Sys.time()
for(k in seq_len(nrow(designs))){
  method <- designs$method[k]
  d <- designs$d[k]
  construction <- designs$construction[k]
  message(
    sprintf(
      "%2d/%d %-8s d = %3d %-10s [%s]",
      k,
      nrow(designs),
      method,
      d,
      construction,
      format(Sys.time() - t0, digits = 3)
    )
  )
  
  results[[length(results) + 1]] <- run_experiment(
    method,
    d,
    construction,
    ns = ns,
    M = M
  )
  
}

results <- do.call(rbind, results)


saveRDS(results, "results/results100.rds")
message("Took ", format(Sys.time() - t0, digits = 3))

################################################################################

t0 <- Sys.time()
effdim <- list()

for(construction in c("sequential", "bb")){
  for(d in dims){
    message(sprintf("effdim: %-10s d = %3d", construction, d))
    effdim[[length(effdim) + 1]] <- unexplained_share(d, construction)
  }
}
effdim <- do.call(rbind, effdim)

saveRDS(effdim, "results/effdim.rds")
message("Effective dimension took ", format(Sys.time() - t0, digits = 3))

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")