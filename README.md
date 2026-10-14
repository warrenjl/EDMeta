# EDMeta: Exposure Distribution Meta-Regression Model

EDMeta implements a hierarchical Bayesian Gaussian process meta-regression framework for synthesizing evidence from environmental epidemiology studies with varying exposure ranges. Each study contributes an estimated effect and standard error together with a summary of its study-specific exposure distribution; the model represents the latent exposure-response function nonparametrically via a Gaussian predictive process indexed over the global exposure axis and integrates that function over each study's exposure distribution to obtain study-specific predictions. The model is fit using Markov chain Monte Carlo sampling techniques. Please see the "EDMeta\_Model\_Details" and "EDMeta\_Example" folders for more specific information regarding the statistical model and package use details, respectively.

## Reference

Warren JL, Xu X, Guo P, and Liew Z (2026). Exposure distribution meta-regression for estimating nonlinear exposure-response relationships. *In Preparation*.

