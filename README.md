# exactnl

`exactnl` is a Stata command that estimates a share-form nested logit demand
model by exact profile likelihood when choice sets vary across markets.

Berry's (1994) linear IV estimator regresses the log share-ratio on the
within-nest log share and instruments it with the count of rival products in
the nest. When nest size is correlated with unobserved mean utility, that count
instrument is invalid and the nesting parameter is biased. `exactnl` instead
writes the likelihood of the observed shares directly, including the
change-of-variables Jacobian of the nested-logit share map. For a nest of
`J` products the Jacobian contributes `(J - 1) * ln(1 - sigma)` to the log
likelihood — the term the linear regression discards. The command concentrates
the linear parameters out by fixed-effect partialling, maximizes the profile
likelihood over the nesting parameters on a grid, and recovers the remaining
coefficients by a final regression at the optimum. Price endogeneity is handled
by a control function (a first-stage `reghdfe` residual added to the structural
equation).

See `help exactnl` after installation for the full syntax, options, and stored
results, and `exactnl_example.do` for worked examples.

## Installation

```stata
net install exactnl, from(https://raw.githubusercontent.com/lisa-george/exactnl/main/)
```

This install line works once the repository is public.

## Verification

`exactnl_cert_public.do` is a self-contained certification script that lets any
user verify the estimator on simulated data. It generates a nested-logit sample
from a fixed seed (no confidential data) and confirms that `exactnl` recovers
the true nesting parameter and matches the paper's simulation-code estimate on
the same grid.

## Citation

If you use this command, please cite:

> Baker, Matthew J. and Lisa M. George. "Demand Estimation with Variable Choice
> Sets: A Likelihood Correction for BLP Models." Working paper.
> SSRN: <SSRN-LINK-PLACEHOLDER>

## Issues

Report problems at https://github.com/lisa-george/exactnl/issues.
