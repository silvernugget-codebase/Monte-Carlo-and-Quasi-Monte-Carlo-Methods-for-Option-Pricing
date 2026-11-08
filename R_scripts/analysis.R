library(tidyverse)

################################################################################

summarise_runs <- function(
  results100,
  payoff = c("arithmetic", "geometric"),
  alpha = 0.05
){
  payoff <- match.arg(payoff)
  z <- qnorm(1 - alpha / 2)
  results100 |> 
    filter(estimator != "QMC") |> 
    mutate(
      value = .data[[payoff]],
      method = factor(
        method,
        levels = c("mc", "halton", "ghalton", "sobol"),
        labels = c("MC", "Halton", "Generalized Halton", "Sobol")
      ),
      construction = factor(
        construction,
        levels = c("sequential", "bb"),
        labels = c("Sequential", "Brownian Bridge")
      )
    ) |> 
    group_by(method, d, construction, n) |> 
    summarise(
      M = length(value),
      est = mean(value),
      sd_n = sd(value), # sd of one n-point estimate across M replicates
      se = sd_n / sqrt(M),# sd of \hat{I}_{n,M}; n*M points
      lwr = est - z * se,
      upr = est + z * se,
      .groups = "drop"
    )
}

################################################################################

# VRF_n = sigma_f^2 / (n * sigma_n^2)

add_vrf <- function(summ){
  sigma2_f <- summ |> 
    filter(method == "MC") |> 
    group_by(d) |> 
    summarise(
      sigma2_f = mean(n * sd_n^2), # variance of a single PO = mean(n * spread
      # of average of  n paths)
      .groups = "drop"
    )
  
  summ |> 
    left_join(sigma2_f, by = "d") |> 
    mutate(
      vrf = sigma2_f / (n * sd_n^2)
    )
}

vrf_table <- function(
  summ_vrf,
  n_fixed = max(summ_vrf$n)
){
  summ_vrf |> 
    filter(n == n_fixed, method != "MC") |> 
    mutate(
      vrf = trimws(formatC(signif(vrf, 2), format = "fg",
                           digits = 2, big.mark = ","))
    ) |> 
    select(
      Construction = construction,
      Point_Set = method,
      d,
      vrf
    ) |> 
    arrange(
      Construction,
      Point_Set,
      d
    ) |> 
    pivot_wider(
      names_from = d,
      values_from = vrf
    )
}

################################################################################

# convergence
#sigma_n = C * n^(- rate) -> log sigma n = log C - rate * log n
# MC: 0.5
# RQMC: close to 1 for low effective dimension

rates <- function(summ){
  summ |> 
    group_by(method, d, construction) |> 
    summarise(
      rate = - unname(coef(lm(log(sd_n) ~ log(n)))[2]),
      .groups = "drop"
    )
}

rates_table <- function(summ){
  rates(summ) |> 
    mutate(rate = sprintf("%.2f", rate)) |> 
    arrange(
      construction,
      method,
      d
    ) |> 
    pivot_wider(
      names_from = d,
      values_from = rate
    )
}

################################################################################

# Coverage of the closed-form price only available for d=1 and geometric
# coverage "by=NULL" = across all methods, dims, construction.
# the finer we split by= the greater the noise

coverage <- function(summ_geo, by = NULL){
  summ_geo |> 
    filter(!(d == 1 & construction == "Brownian Bridge")) |>
    # avoid a double reading since d=1 + bb == d=1 + seq
    mutate(
      exact = sapply(d, geometric_price), # applies geometric_price on all d
      covered = lwr <= exact & exact <= upr # T / F
    ) |> 
    group_by(across(all_of(by))) |> 
    summarise(
      share_covered = mean(covered), # frac of T
      intervals = length(covered), # 
      .groups = "drop"
    )
}

################################################################################

# price table: arithmetic asian call at N = n*M total points, 
# arithmetic price must be greater than geometric

price_table <- function(
  summ_arith,
  n_fixed = max(summ_arith$n),
  rqmc = "Sobol",
  cons = "Brownian Bridge"
){
  f4 <- function(x)
    formatC(x, format = "f", digits = 6)
  ci <- function(est, se)
    paste0("$", f4(est), " \\pm", f4(qnorm(0.975) * se), "$")
  # estimate +- z_(1-0.05/2) * se
  
  mc_rows <- summ_arith |> 
    filter(
      n == n_fixed,
      method == "MC",
      construction == "Sequential" # identical to bridge
    )
  qmc_rows <- summ_arith |> 
    filter(
      n == n_fixed,
      method == rqmc,
      construction == cons
    )
  
  tibble(
    d = mc_rows$d,
    geometric = paste0("$", f4(sapply(mc_rows$d, geometric_price)), "$"),
    mc = ci(mc_rows$est, mc_rows$se),
    rqmc = ci(qmc_rows$est, qmc_rows$se)
  )
    
}

################################################################################

# Validation plot

plot_validation <- function(
  summ_geo,
  dim,
  cons = "Brownian Bridge"
){
  summ_geo |> 
    filter(d == dim, construction == cons) |> 
    ggplot(aes(n, est, color = method, fill = method)) +
    geom_ribbon(
      aes(ymin = lwr, ymax = upr),
      alpha = 0.25,
      color = NA
    ) +
    geom_line() +
    geom_hline(
      yintercept = geometric_price(dim),
      linetype = "dotted"
    ) +
    scale_x_continuous(
      trans = "log2",
      breaks = 2^c(8, 10, 12, 14, 16),
      labels = c("$2^8$", "$2^{10}$", "$2^{12}$", "$2^{14}$", "$2^{16}$")
    ) +
    labs(
      x = "$n$",
      y = "Price Estimate",
      color = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom")
}

################################################################################

# Error plot
plot_error <- function(summ, dim){
  df <- summ |> 
    filter(d == dim)
  n0 <- min(df$n)
  
  ref <- df |> 
    group_by(construction) |> 
    mutate(
      ref_mc = sd_n[method == "MC" & n == n0] * sqrt(n0 / n)
    ) |> 
    ungroup() |> 
    filter(method == "MC")
  
  ggplot(df, aes(n, sd_n, color = method)) +
    geom_line(
      data = ref,
      aes(n, ref_mc),
      linetype = "dashed",
      color = "grey40",
      inherit.aes = FALSE
    ) +
    geom_line(linewidth = 0.3) +
    geom_point(size = 0.8) +
    facet_wrap(~ construction) +
    scale_x_continuous(
      trans = "log2",
      breaks = 2^c(8, 10, 12, 14, 16),
      labels = c("$2^8$", "$2^{10}$", "$2^{12}$", "$2^{14}$", "$2^{16}$")
    ) +
    scale_y_log10() +
    labs(
      x = "$n$",
      y = "Standard Deviation of One Estimate",
      color = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "bottom",
      panel.spacing.x = unit(1.5, "lines")
    )
}

# VRF plot
plot_vrf <- function(
  summ_vrf,
  n_fixed = max(summ_vrf$n)
){
  summ_vrf |> 
    filter(n == n_fixed, method != "MC") |> 
    ggplot(aes(d, vrf, color = method)) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    geom_line() +
    geom_point(size = 1) +
    facet_wrap(~ construction) +
    scale_x_log10(breaks = c(1, 4, 12, 52, 252)) +
    scale_y_log10() +
    labs(
      x = "$d$",
      y = "VRF",
      color = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "bottom",
      panel.spacing.x = unit(1.5, "lines")
    )
}

################################################################################

# effective dimension

truncation_level <- 0.99

truncation_dim <- function(ed, level = truncation_level){
  ed |> 
    group_by(d, construction) |> 
    summarise(
      d_t = min(s[(1 - share) >= level]),
      .groups = "drop"
    )
}

plot_effective_dim <- function(ed, level = truncation_level){
  ed |> 
    filter(share > 0) |> 
    ggplot(aes(s, share, color = construction)) +
    geom_hline(
      yintercept = 1 - level,
      linetype = "dashed"
    ) +
    geom_line() +
    facet_wrap(
      ~ d,
      labeller = as_labeller(function(x) paste0("$d =", x, "$"))
    ) +
    scale_x_log10() +
    scale_y_log10(
      labels = function(x)
        ifelse(is.na(x), NA, sprintf("$10^{%d}$", as.integer(round(log10(x)))))
    ) +
    labs(
      x = "$s$ (Number of Leading Coordinates)",
      y = "Unexplained Variance Share",
      color = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom")
}



