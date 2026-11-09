# Monte Carlo vs. Quasi-Monte Carlo for Asian Option Pricing

This repository contains the code for my Bachelor's thesis titled *Monte Carlo and Quasi-Monte Carlo Integration: Pricing Discretely Monitored Arithmetic Asian Call Options in the Black-Scholes Model*. It contains a file that exactly reproduces the numerical experiment in the final section of the thesis and the resulting figures and tables.

The research question asks

> *For a discretely monitored arithmetic Asian call on a non-dividend-paying stock in the Black-Scholes model, how does (Randomized) Quasi-Monte Carlo compare with standard Monte Carlo as the dimension of the integral grows and how much of the difference depends on how the Brownian path is constructed (sequential vs. Brownian bridge)?*

The thesis compares the standard Monte Carlo (MC) method with Randomized Quasi-Monte Carlo (RQMC) for pricing an arithmetic Asian call with up to 252 monitoring dates for two types of Brownian path construction: The random walk construction and the Brownian bridge construction. The closed-form solution of the geometric Asian call option served as the validation benchmark to ensure a correct implementation. Using the Brownian bridge construction, RQMC drastically reduces the variance over MC; under the random walk construction, this advantage shrinks sharply as the number of monitoring dates increases.

**Reproducing the results:**

1. Install and load the following `R` packages:
   - `qrng`: the point sets (Halton, generalized Halton and Sobol).
   - `withr`: freezes the random digital shift of `qrng::ghalton()` with a fixed seed.
   - `tidyverse`: data wrangling and plots in **`analysis.R`**.
   - `knitr`: runs the `R` code when rendering.
   - `kableExtra`: formats the tables.
   - `tikzDevice`: renders the figures in LaTeX.
   
2. Run **`R_scripts/experiments.R`** from the project root. It calls **`run_experiment()`** for every combination of point set, dimension and construction and saves the results to `results/results100.rds`, then calls **`unexplained_share()`** for every dimension and both constructions and saves the unexplained shares to `results/effdim.rds` (about 105 minutes in total).

3. Render `Results.qmd` with Quarto, which sources **`functions.R`** and **`analysis.R`**. This requires a LaTeX installation, since the PDF is built with `pdflatex` and the figures with `tikzDevice`.

**Important Note:** The `results/` folder already contains `results100.rds` and `effdim.rds`, hence `Results.qmd` can be rendered (step 3.) without the need to re-run **`run_experiment()`**, provided all packages from step 1. are loaded.

All random numbers are seeded to make the results exactly reproducible with the same software versions. The code was writtin under R version 4.6.0 and `qrng` version 0.0-11.


