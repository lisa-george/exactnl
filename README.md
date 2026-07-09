# exactnl

exactnl is a Stata command that estimates the share-form nested logit
demand model by exact maximum likelihood.

## What it does

The share-form nested logit is estimated on a transformed dependent
variable: observed shares are mapped to the log share-ratio
ln(s_j) - ln(s_0) following Berry (1994). A likelihood written for the
transformed variable must include the Jacobian of that transformation, and
in the nested logit the Jacobian depends on the nesting parameter and the
nest size: each nest of J products contributes (J - 1) * ln(1 - sigma) to
the log likelihood. Standard practice estimates the share equation as a
linear regression and omits this term. When every market has the same
number of products the term is constant and its omission is harmless. When
choice sets vary across markets, omission biases the nesting parameter,
and the common remedy, instrumenting the within-nest share with product
counts, fails whenever product counts also enter demand directly.

exactnl maximizes the exact likelihood instead. It concentrates the linear
parameters out by fixed-effect partialling, profiles the likelihood over
the nesting parameters, and recovers the remaining coefficients by a final
regression at the optimum. Price endogeneity is handled by a control
function. No instrument for the within-nest share is required, and a
product-count term in mean utility can be estimated directly.

See `help exactnl` after installation for syntax, options, and stored
results, and `exactnl_example.do` for worked examples.

## Installation

    net install exactnl, from(https://raw.githubusercontent.com/lisa-george/exactnl/main/)

## Verification

`exactnl_cert_public.do` is a self-contained certification script. It
simulates nested logit data from a fixed seed and confirms that exactnl
recovers the true nesting parameter. No confidential data required.

## Citation

If you use this command, please cite:

Baker, Matthew J. and Lisa M. George. "Demand Estimation with Variable
Choice Sets: A Likelihood Correction for Nested Logit." Working paper.
[SSRN link]

## Issues

Report problems at https://github.com/lisa-george/exactnl/issues.
