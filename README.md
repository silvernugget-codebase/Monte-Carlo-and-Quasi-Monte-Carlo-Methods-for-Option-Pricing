# Monte Carlo and Quasi-Monte Carlo Integration: Pricing Discretely Monitored Arithmetic Asian Call Options in the Black--Scholes Model

This repository contains the codebase for my Bachelor's thesis titled *Monte Carlo and Quasi-Monte Carlo Integration: Pricing Discretely monitored Arithmetic Asian Call Options in the Black--Scholes Model*. It also contains a file that exactly reproduces the numerical experiment ran at the final section of the thesis and the resulting figures and tables.

The thesis itself is **not** yet included and will be added after the grading is done.

The research question asks

> *For a discretely monitored arithmetic Asian call on a non-dividend-paying stock in the Black--Scholes model, how does (Randomized) Quasi-Monte Carlo compare with standard Monte Carlo as the dimension of the integral grows and how much of the difference depends on how the Brownian path is constructed (sequential vs. Brownian bridge)?*

# Code 

Everything is written in **R**. The code consists of three files in the folder `R_scripts/`.

**`functions.R`**

The building blocks of the estimators:

- **`base_params`**: the base-case (and default) parameters $S_0, K, r, \sigma, T$.
- **`check_model()`**: checks the model parameters and stops with an error message if they are invalid.
- **`qmc_points()`**: generates the first $n$ points of the $d$-dimensional Halton, generalized Halton or Sobol sequence.
- **`ghalton_points()`**: wraps `qrng::ghalton()` for both Halton variants and freezes its built-in random digital shift with a fixed seed.
- **`sobol_points()`**: wraps `qrng::sobol()`; unscrambled by default, so the first point is the origin.
- **`cp_rotation()`**: applies the random shift modulo 1 and moves values that land exactly on $0$ to $2^{-53}$, so that `qnorm()` stays finite.
- **`paths()`**: turns standard normal variates into Brownian paths, sequentially and via the Brownian bridge.
- **`payoff()`**: computes the discounted arithmetic and geometric Asian call payoffs of each path.
- **`geometric_price()`**: computes the closed-form price of the geometric Asian call, the validation benchmark.
- **`unexplained_share()`**: estimates the share of the payoff variance not explained by the first $s$ coordinates.


**`experiments.R`**

Sets the experimental design, (`dims`, `ns`, `M`), defines one function and runs all experiments (see *Reproducing the results* below):

- **`run_experiment()`**: for one point set, dimension and path construction, generates $M$ independent replicates (independent draws for MC, independent random shifts of one fixed point set for RQMC) and returns the arithmetic and geometric price estimates at every $n = 2^8,\dots,2^{16}$.


**`analysis.R`**

Summarizes the runs and produces all tables and figures:

- **`summarise_runs()`**: aggregates the replicates into the estimate, the standard deviation of an $n$-point estimate, the standard error and the 95% confidence interval for each point set, dimension, construction and $n$.
- **`add_vrf()`**: adds the variance reduction factor, with $\sigma_f^2$ estimated from the MC runs.
- **`vrf_table()`**: tabulates the VRFs at $n = 2^{16}$.
- **`rates()`** and **`rates_table()`**: fit and tabulate the empirical convergence rates (negative log-log slope of $\hat{\sigma}_n$ against $n$).
- **`coverage()`**: computes the share of confidence intervals that contain the exact geometric price.
- **`price_table()`**: tabulates the exact geometric price and the arithmetic MC and RQMC prices at $n = 2^{16}$.
- **`plot_validation()`**: plots the estimates with their confidence bands against the exact price.
- **`plot_error()`**: plots $\hat{\sigma}_n$ against $n$.
- **`plot_vrf()`**: plots the VRFs against the dimension.
- **`truncation_dim()`**: computes the truncation dimension $d_t$ from the unexplained shares (99% level).
- **`plot_effective_dim()`**: plots the unexplained share against the number of leading coordinates.


**Reproducing the results:**

1. Install and load the following `R` packages:
   - `qrng`: the point sets (Halton, generalized Halton and Sobol).
   - `withr`: freezes the random digital shift of `qrng::ghalton()` with a fixed seed.
   - `tidyverse`: data wrangling and plots in **`analysis.R`**.
   - `knitr`: runs the `R` code when rendering.
   - `kableExtra`: formats the tables.
   - `tikzDevice`: renders the figures in Latex.
   
2. Run **`R_scripts/experiments.R`** from the project root. It calls **`run_experiment()`** for every combination of point set, dimension and construction and saves the results to `results/results100.rds`, then calls **`unexplained_share()`** for every dimension and both constructions and saves the unexplained shares to `results/effdim.rds` (about 105 minutes in total).

3. Render `Results.qmd` with **Quarto**, which sources **`functions.R`** and **`analysis.R`**. This requires a Latex installation, since the PDF is built with `pdflatex` and the figures with `tikzDevice`.

**Important Note:** The `results/` folder already contains `results100.rds` and `effdim.rds`, hence `Results.qmd` can be rendered (step 3.) without the need to re-run **`run_experiment()`**, provided all packages from step 1. are loaded.

All random numbers are seeded to make the results exactly reproducible with the same software versions. We used R version 4.6.0 and `qrng` version 0.0-11.


