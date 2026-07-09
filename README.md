# exactnl

exactnl is a Stata command that estimates the share-form nested logit
demand model by exact maximum likelihood.

## What it does

The share-form nested logit is estimated on a transformed dependent
variable: observed shares are mapped to the log share-ratio
ln(s_j) - ln(s_0) following Berry (1994). A likelihood written for the
transformed variable must include the Jacobian of that transformation.  In
the nested logit, the Jacobian depends on the choice set and the nesting
parameter.  Each  nest of J products contributes (J - 1) * ln(1 - sigma) to
the log likelihood.  Applied work omits this term. The omission has no effect
when choice sets are fixed across markets.  But when choice sets vary, omitting the Jacobian 
biases estimates of within-group substitution, which feed directly into
elasticities, markups, and welfare calculations.

Guides for structural demand estimation recommend instrumenting the
within-nest share with counts of rival products. Under the exact
likelihood, count instruments are not valid because the same product-count variation that makes them
strong instruments enters the estimating objective directly. However with the corrected likelihood, instruments for the within-nest share are not needed because the dependence between the within-nest share and the structural
error is exactly what the Jacobian accounts for.

exactnl maximizes the exact likelihood. It concentrates the linear
parameters out by fixed-effect partialling, profiles the likelihood over
the nesting parameters, and recovers the remaining coefficients at the
optimum. Price endogeneity is handled by a control function. A
product-count term in mean utility can be estimated directly, so crowding
of the product space (Ackerberg and Rysman 2005) is estimable rather than
assumed away.

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
