
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`medshift`

[![Travis-CI Build
Status](https://travis-ci.org/nhejazi/medshift.svg?branch=master)](https://travis-ci.org/nhejazi/medshift)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/medshift?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/medshift)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/medshift/master.svg)](https://codecov.io/github/nhejazi/medshift?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Causal Mediation Analysis for Stochastic Interventions

**Authors:** [Nima Hejazi](https://nimahejazi.org) and [Iván
Díaz](https://idiaz.xyz)

-----

## What’s `medshift`?

The `medshift` R package is designed to provide facilities for
estimating a parameter that arises in a decomposition of the population
intervention causal effect into the (in)direct effects under stochastic
interventions in the setting of mediation analysis. `medshift` is
designed as an implementation to accompany the methodology described in
Dı́az and Hejazi (2020). Implemented estimators include the classical
substitution (G-computation) estimator, an inverse probability weighted
(IPW) estimator, an efficient one-step estimator using cross-fitting
(Pfanzagl and Wefelmeyer 1985; Zheng and van der Laan 2011; Chernozhukov
et al. 2018), and a cross-validated targeted minimum loss (TML)
estimator (van der Laan and Rose 2011; Zheng and van der Laan 2011).
`medshift` integrates with the [`sl3` R
package](https://github.com/tlverse/sl3) (Coyle et al. 2022) to allow
constructed estimators to leverage machine learning for nuisance
estimation.

-----

## Installation

Install the *most recent version* from the `master` branch on GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("nhejazi/medshift")
```

-----

## Example

To illustrate how `medshift` may be used to estimate the effect of
applying a stochastic intervention to the treatment (`A`) while keeping
the mediator(s) (`Z`) fixed, consider the following example:

``` r
library(data.table)
library(medshift)

# produces a simple data set based on ca causal model with mediation
make_simple_mediation_data <- function(n_obs = 1000) {
  # baseline covariate -- simple, binary
  W <- rbinom(n_obs, 1, prob = 0.50)

  # create treatment based on baseline W
  A <- as.numeric(rbinom(n_obs, 1, prob = W / 4 + 0.1))

  # single mediator to affect the outcome
  z1_prob <- 1 - plogis((A^2 + W) / (A + W^3 + 0.5))
  Z <- rbinom(n_obs, 1, prob = z1_prob)

  # create outcome as a linear function of A, W + white noise
  Y <- Z + A - 0.1 * W + rnorm(n_obs, mean = 0, sd = 0.25)

  # full data structure
  data <- as.data.table(cbind(Y, Z, A, W))
  setnames(data, c("Y", "Z", "A", "W"))
  return(data)
}

# set seed and simulate example data
set.seed(75681)
example_data <- make_simple_mediation_data()

# compute one-step estimate for an incremental propensity score intervention
# that triples (delta = 3) the individual-specific odds of receiving treatment
os_medshift <- medshift(W = example_data$W, A = example_data$A,
                        Z = example_data$Z, Y = example_data$Y,
                        delta = 3, estimator = "onestep",
                        estimator_args = list(cv_folds = 3))
summary(os_medshift)
#>      lwr_ci   param_est      upr_ci   param_var    eif_mean   estimator 
#>      0.7401    0.788136    0.836172    0.000601 1.64686e-17     onestep
```

For details on how to use data adaptive regression (machine learning)
techniques in the estimation of nuisance parameters, consider consulting
the vignette that accompanies this package.

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/medshift/issues).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/nhejazi/medshift/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `medshift` R package, please cite the following:

``` 
    @article{diaz2020causal,
      title={Causal mediation analysis for stochastic interventions},
      author={D{\'\i}az, Iv{\'a}n and Hejazi, Nima S},
      year={2020},
      url = {https://doi.org/10.1111/rssb.12362},
      doi = {10.1111/rssb.12362},
      journal={Journal of the Royal Statistical Society: Series B
        (Statistical Methodology)},
      volume={},
      number={},
      pages={},
      publisher={Wiley Online Library}
    }

    @manual{hejazi2020medshift,
      author = {Hejazi, Nima S and D{\'\i}az, Iv{\'a}n},
      title = {{medshift}: Causal mediation analysis for stochastic
        interventions},
      year  = {2020},
      url = {https://github.com/nhejazi/medshift},
      note = {R package version 0.1.4}
    }
```

-----

## License

© 2018-2022 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2018-2022 Nima S. Hejazi
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

-----

## References

<div id="refs" class="references">

<div id="ref-chernozhukov2018double">

Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo,
Christian Hansen, Whitney Newey, and James Robins. 2018.
“Double/Debiased Machine Learning for Treatment and Structural
Parameters.” *The Econometrics Journal* 21 (1).
<https://doi.org/10.1111/ectj.12097>.

</div>

<div id="ref-coyle-gh-sl3">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, Rachael V Phillips, and
Oleg Sofrygin. 2022. *sl3: Modern Pipelines for Machine Learning and
Super Learning*. <https://github.com/tlverse/sl3>.
<https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-diaz2020causal">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)*.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-pfanzagl1985contributions">

Pfanzagl, J, and W Wefelmeyer. 1985. “Contributions to a General
Asymptotic Statistical Theory.” *Statistics & Risk Modeling* 3 (3-4):
379–88.

</div>

<div id="ref-vdl2011targeted">

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

</div>

<div id="ref-zheng2011cross">

Zheng, Wenjing, and Mark J van der Laan. 2011. “Cross-Validated Targeted
Minimum-Loss-Based Estimation.” In *Targeted Learning*, 459–74.
Springer.

</div>

</div>
