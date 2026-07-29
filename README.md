# exactnl

exactnl is a Stata command that estimates the share-form nested logit
demand model by maximum likelihood, using the exact density of the
transformed dependent variable.

## What it does

The share-form nested logit is estimated on a transformed dependent
variable: observed shares are mapped to the log share-ratio
ln(s_j) - ln(s_0) following Berry (1994). A likelihood written for the
transformed variable must include the Jacobian of that transformation.  In
the nested logit, the Jacobian depends on the choice set and the nesting
parameter.  Each nest of J products contributes (J - 1) * ln(1 - sigma) to
the log likelihood.  Applied work omits this term. Because the term depends
on the nesting parameter, omitting it biases likelihood estimates of
within-group substitution, which feed directly into elasticities, markups,
and welfare calculations. The bias is entangled with product counts when
choice sets vary across markets.

Guides for structural demand estimation recommend instrumenting the
within-nest share with counts of rival products. Under the corrected
likelihood such instruments are not needed: the dependence between the
within-nest share and the structural error is exactly what the Jacobian
accounts for. And when product counts enter demand directly, as when
products crowd the product space, count instruments fail exclusion and are
not valid. The corrected likelihood makes count instruments unnecessary;
crowding makes them invalid.

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

## Requirements

Stata 17 or later. Requires `reghdfe` and `ftools`:

    ssc install reghdfe
    ssc install ftools

## Syntax

    exactnl depvar [indepvars] [(endogvar = instvars)], nest(varname) market(varlist) [options]

`depvar` is the user-constructed dependent variable ln(s_j) - ln(s_0).
Within-nest shares are supplied either as raw shares via `share()` and
`outside()`, or precomputed via `wnshare()`.

Two nests, county and year fixed effects, price instrumented by a cost
shifter through a control function, a common crowding term, and clustered
standard errors:

    exactnl y x1 x2 (price = costshifter), nest(nestid) market(county year) ///
        absorb(county year) share(s_j) outside(s_0) crowding("lnJ") cluster(county)

See `help exactnl` for the full syntax, options, and stored results.

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
