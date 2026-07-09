# exactnl

exactnl is a Stata command that estimates the share-form nested logit demand
model by exact maximum likelihood.

## What it does

The Berry (1994) inversion turns nested logit demand into a linear equation:
the log share-ratio regressed on characteristics, price, and the within-nest
log share. The exact likelihood of this equation contains a Jacobian term
from the share transformation: each nest of J products contributes
(J - 1) * ln(1 - sigma) to the log likelihood. Standard practice estimates
the linear equation by IV and discards this term. exactnl includes it. The
command concentrates the linear parameters out by fixed-effect partialling,
maximizes the profile likelihood over the nesting parameters, and recovers
the remaining coefficients by a final regression at the optimum. Price
endogeneity is handled by a control function. No instrument for the
within-nest share is required.

## When the correction matters

If every market has the same number of products, the Jacobian term is a
constant and omitting it is harmless: exactnl and the standard estimator
agree.

If choice sets vary across markets, the term varies with the nesting
parameter and its omission biases the estimate. The common remedy in
applied work is to instrument the within-nest share with the count of
products in the nest. This works when product counts affect only how shares
divide within a nest. It fails when product counts also enter demand
directly, as when products crowd a congested product space (Ackerberg and
Rysman 2005): the count then belongs in the demand equation and cannot be
excluded. exactnl requires no share instrument, so it remains valid in
exactly these settings, and it permits a product-count term in mean utility
so that crowding can be estimated rather than assumed away.

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

Baker, Matthew J. and Lisa M. George. "Maximum Likelihood Estimation of the
Share-Form Nested Multinomial Logit Model." Working paper. [SSRN link]

## Issues

Report problems at https://github.com/lisa-george/exactnl/issues.
