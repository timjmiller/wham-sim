## WHAM description + simulation tests

*The Woods Hole Assessment Model (WHAM): a general state-space assessment framework that incorporates time- and age-varying processes via random effects and links to environmental covariates*

Brian C. Stock and Timothy J. Miller
Population Dynamics Branch, NEFSC, NOAA Fisheries

### Abstract

The rapid changes observed in many marine ecosystems that support fisheries pose a challenge to stock assessment and management predicated on time-invariant productivity and considering species in isolation. In single-species assessments, two main approaches have been used to account for productivity changes: allowing biological parameters to vary stochastically over time (empirical), or explicitly linking population processes such as recruitment (*R*) or natural mortality (*M*) to environmental covariates (mechanistic). Here, we describe the Woods Hole Assessment Model (WHAM) framework and software package, which combines these two approaches. WHAM can estimate time- and age-varying random effects on annual transitions in numbers at age (NAA), *M*, and selectivity, as well as fit environmental time-series with process and observation errors, missing data, and nonlinear links to *R* and *M*. WHAM can also be configured as a traditional statistical catch-at-age (SCAA) model in order to easily bridge from status quo models and test them against models with state-space and environmental effects, all within a single framework.

We fit models with and without (independent or autocorrelated) random effects on NAA, *M*, and selectivity to data from five stocks with a broad range of life history, fishing pressure, number of ages, and time-series length. Models that included random effects performed well across stocks and processes, especially random effects models with a two dimensional (2D) first-order autoregressive (AR1) covariance structure over age and year. We conducted simulation tests and found negligible or no bias in estimation of important assessment outputs (SSB, *F*, stock status, and catch) when the operating and estimation models matched. However, bias in SSB and *F* was often non-trivial when the estimation model was less complex than the operating model, especially when models without random effects were fit to data simulated from models with random effects. Bias of the variance and correlation parameters controlling random effects was also negligible or slightly negative as expected. Our results suggest that WHAM can be a useful tool for stock assessment when environmental effects on *R* or *M*, or stochastic variation in NAA transitions, *M*, or selectivity are of interest. In the U.S. Northeast, where the productivity of several groundfish stocks has declined, conducting assessments in WHAM with time-varying processes via random effects or environment-productivity links may account for these trends and potentially reduce retrospective bias.

### Keywords

state-space; stock assessment; random effects; time-varying; environmental effects; recruitment; natural mortality; Template Model Builder (TMB)
