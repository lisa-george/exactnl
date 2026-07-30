{smcl}
{* *! version 1.0.0  29jul2026  Baker & George}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "nlogit" "help nlogit"}{...}
{viewerjumpto "Syntax" "exactnl##syntax"}{...}
{viewerjumpto "Description" "exactnl##description"}{...}
{viewerjumpto "Do not instrument the within-nest share" "exactnl##warning"}{...}
{viewerjumpto "Options" "exactnl##options"}{...}
{viewerjumpto "Stored results" "exactnl##results"}{...}
{viewerjumpto "Examples" "exactnl##examples"}{...}
{viewerjumpto "Remarks" "exactnl##remarks"}{...}
{viewerjumpto "Authors" "exactnl##authors"}{...}
{viewerjumpto "References" "exactnl##references"}{...}
{title:Title}

{phang}
{bf:exactnl} {hline 2} Share-form nested logit by maximum likelihood with the
exact change-of-variables correction, fixed-effect absorption, and a
control-function price instrument

{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:exactnl} {it:depvar} [{it:indepvars}] [{cmd:(}{it:endogvar} {cmd:=} {it:instvars}{cmd:)]}{cmd:,}
{opth nest(varname)} {opth market(varlist)}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model (required)}
{synopt:{opth nest(varname)}}categorical variable defining the nests (1..G per
market); one nesting parameter sigma is estimated per {it:level} of {it:nest()}{p_end}
{synopt:{opth market(varlist)}}market identifiers (e.g. {cmd:county year}); a
nest is a {it:market x nest()} cell{p_end}

{syntab:Within-nest shares (supply ONE of the two)}
{synopt:{opth share(varname)}}product-level market share s_j; combined with
{cmd:outside()} the command builds ln s_(j|g) and validates the shares{p_end}
{synopt:{opth out:side(varname)}}outside-good share s_0 (required with
{cmd:share()}){p_end}
{synopt:{opth wns:hare(varname)}}precomputed within-nest log share ln s_(j|g),
for data that do not carry raw shares{p_end}

{syntab:Specification}
{synopt:{opth absorb(varlist)}}fixed effects to absorb, reghdfe-style (e.g.
{cmd:county year}){p_end}
{synopt:{opt crow:ding(str)}}{cmd:lnJ} adds tau*ln(J_gm) to mean utility with a
common tau; {cmd:lnJ, bynest} adds a nest-specific tau_g{p_end}
{synopt:{opt noJAC:obian}}omit the Jacobian term (for comparison; prints a
warning that sigma is biased upward){p_end}

{syntab:Inference}
{synopt:{opth cluster(varname)}}cluster-robust (sandwich) SEs for the linear
parameters{p_end}
{synopt:{opt boot:strap(reps)}}cluster block-bootstrap SEs for sigma_g, tau and
the linear parameters; use with {cmd:cluster()}{p_end}

{syntab:Numerical}
{synopt:{opt grid(lo hi step)}}profile grid for sigma (default {cmd:0.001 0.999
0.001}){p_end}
{synoptline}

{p 4 6 2}
{it:depvar} must be the user-constructed logit dependent variable
ln(s_j) - ln(s_0).  The parenthesised {cmd:(}{it:endogvar} {cmd:=} {it:instvars}{cmd:)}
group requests a control-function first stage for the (endogenous) price; omit it
to treat all regressors as exogenous.

{p 4 6 2}
{cmd:exactnl} requires Stata 17 or later and {help reghdfe} (which in turn
requires {cmd:ftools}).  Install both with {cmd:ssc install reghdfe} and
{cmd:ssc install ftools}.

{marker description}{...}
{title:Description}

{pstd}
{cmd:exactnl} estimates the share-form nested logit demand model by exact
maximum likelihood.  The model is estimated on a transformed dependent
variable: observed shares are mapped to the log share-ratio
ln(s_j) - ln(s_0) following Berry (1994).  A likelihood written for the
transformed variable must include the Jacobian of that transformation, and in
the nested logit the Jacobian depends on the choice set and the nesting
parameter: each nest of J_gm products contributes
{bf:(J_gm - 1) ln(1 - sigma_g)} to the log likelihood.  Applied work omits this
term.  Because the term depends on sigma_g, its omission biases likelihood
estimates of within-group substitution; when choice sets vary across markets
the bias is entangled with product counts.  {cmd:exactnl} includes the term.
Concentrating the linear parameters out by fixed-effect partialling, the
profile log likelihood is

{p 8 8 2}
LL(sigma) = -N/2 ln(SSR(sigma)/N) + sum_g [ sum_m (J_gm - 1) ] ln(1 - sigma_g),

{pstd}
where SSR(sigma) is a quadratic form in sigma built once from the FE-partialled
residuals of {it:depvar} and of the within-nest log shares interacted with each
nest.  {cmd:exactnl} maximises LL over sigma in [lo,hi]^G (an exhaustive grid
for G<=2; coordinate golden-section ascent for G>=3), then recovers the linear
parameters (alpha on price, beta on {it:indepvars}, and tau on the crowding
term) by one final FE regression at the optimum.

{pstd}
No instrument for the within-nest share is used.  Applied work commonly
instruments ln s_(j|g) with counts of rival products in the nest; under the
exact likelihood such instruments are not needed.  The dependence between the
within-nest share and the structural error is exactly what the Jacobian
accounts for.  When product counts also enter mean utility directly, as when
products crowd the product space, count instruments fail the exclusion
restriction and are not valid.  The corrected likelihood makes count
instruments unnecessary; crowding makes them invalid.  Because no share
instrument is required, a product-count (crowding) term in mean utility can be
estimated directly via {cmd:crowding()}, in the spirit of
Ackerberg and Rysman (2005).

{pstd}
Price endogeneity is a separate matter and still requires instruments.  It is
handled by a control function: the first stage regresses {it:endogvar} on the
excluded instrument(s), the exogenous regressors, the crowding term and the
absorbed fixed effects (via {help reghdfe}), and the first-stage residual is
added to the structural equation.  Standard errors on the linear parameters
are the {cmd:cluster()} sandwich from the final regression; the sigma_g
standard errors are profile-curvature (observed-information) approximations,
and a cluster block bootstrap ({cmd:bootstrap()}) is recommended for inference
on sigma_g.  The command also reports the boundary likelihood-ratio test of
H0: all sigma_g = 0 against the one-sided alternative.  The printed 5%
critical value is the Kodde-Palm (1986) least-favourable upper bound (2.71 for
G=1, 5.14 for G=2), which is valid regardless of the dependence among the
sigma_g estimates; the exact chi-bar-squared mixture quantile (Self and Liang 1987; valid under
independence of the sigma_g estimates) is stored in {cmd:e(LR_crit5_exact)}
(see {help exactnl##results:Stored results}).

{pstd}
The Jacobian term is implied by the change of variables and holds for any
absolutely continuous error density.  The implementation is Gaussian, so
{cmd:exactnl} is exact maximum likelihood under normal structural errors and a
quasi-maximum-likelihood estimator otherwise.  The point estimate of sigma_g is
robust to heavy-tailed errors in simulation; the curvature-based standard
errors are not.  Use {cmd:cluster()} and {cmd:bootstrap()} for inference.

{marker warning}{...}
{title:Do not instrument the within-nest share}

{pstd}
{bf:The Jacobian already accounts for the dependence between the within-nest
share and the structural error.}  Supplying a within-share instrument or
control function in addition to the correction applies the same adjustment
twice and biases sigma_g downward.  The parenthesised
{cmd:(}{it:endogvar} {cmd:=} {it:instvars}{cmd:)} group is for price only.
Do not place ln s_(j|g) on the left of the equals sign, and do not include
rival product counts among {it:instvars}.

{pstd}
Product counts may appear in mean utility through {cmd:crowding()}.  That is a
regressor, not an instrument, and is the intended way to let nest size affect
demand.

{marker options}{...}
{title:Options}

{phang}
{opth nest(varname)} (required) is the categorical nest variable.  Its distinct
levels index the nesting parameters sigma_1..sigma_G (sorted ascending); a nest
is a {it:market() x nest()} cell.

{phang}
{opth market(varlist)} (required) identifies markets.  J_gm and the within-nest
shares are computed within {it:market() x nest()} cells.

{phang}
{opth share(varname)} / {opth outside(varname)} supply the raw product share s_j
and outside share s_0.  The command builds ln s_(j|g) = ln(s_j) - ln(sum of s_j
in the nest) and validates that shares lie in (0,1) and that within-market share
totals do not exceed 1.

{phang}
{opth wnshare(varname)} supplies a precomputed within-nest log share ln s_(j|g)
directly, for data that store the within-nest log share but not the raw shares.
Exactly one of {cmd:wnshare()} or {cmd:share()}+{cmd:outside()} must be given.

{phang}
{opth absorb(varlist)} lists the fixed effects to absorb ({help reghdfe} syntax).

{phang}
{opt crowding(str)} adds a nest-size (crowding) term to mean utility.
{cmd:crowding(lnJ)} adds tau*ln(J_gm) with a single common tau;
{cmd:crowding(lnJ, bynest)} adds a separate tau_g for each nest.

{phang}
{opt noJACobian} drops the (J_gm-1)ln(1-sigma_g) term, reducing the objective to
the Gaussian sum of squares.  The resulting sigma_g are biased upward; a
warning is printed.

{phang}
{opth cluster(varname)} requests cluster-robust SEs for the linear parameters.

{phang}
{opt bootstrap(reps)} runs a cluster block bootstrap of the full estimator,
resampling clusters, rebuilding the nest counts on the bootstrap replica id, and
refitting; it reports bootstrap SEs for sigma_g, tau and the linear parameters.
Requires {cmd:cluster()} to be the leading {cmd:market()} variable.

{phang}
{opt grid(lo hi step)} sets the profile grid (default {cmd:0.001 0.999 0.001}).

{marker results}{...}
{title:Stored results}

{pstd}{cmd:exactnl} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_markets)}}number of markets{p_end}
{synopt:{cmd:e(G)}}number of nests (distinct levels of {cmd:nest()}){p_end}
{synopt:{cmd:e(LL)}}corrected profile log likelihood at the optimum{p_end}
{synopt:{cmd:e(LL_sigma0)}}log likelihood at sigma = 0{p_end}
{synopt:{cmd:e(LR)}}boundary LR statistic 2[e(LL)-e(LL_sigma0)]{p_end}
{synopt:{cmd:e(LR_crit5)}}Kodde-Palm (1986) 5% least-favourable upper bound for
the boundary LR test, c solving 0.5 chi2tail(G-1,c) + 0.5 chi2tail(G,c) = 0.05
(2.71 for G=1; 5.14 for G=2); conservative and valid regardless of the
dependence among the sigma_g estimates -- this is the value printed{p_end}
{synopt:{cmd:e(LR_crit5_exact)}}exact 5% quantile of the chi-bar-squared mixture
sum_k C(G,k)2^(-G) chi2_k (2.71 for G=1; ~4.23 for G=2); exact only when the
sigma_g estimates are independent{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:exactnl}{p_end}
{synopt:{cmd:e(depvar)}}name of {it:depvar}{p_end}
{synopt:{cmd:e(crowding)}}the {cmd:crowding()} specification, if any{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector: {it:indepvars} beta, alpha (price),
tau/tau_g, then sigma1..sigmaG{p_end}
{synopt:{cmd:e(V)}}variance matrix (cluster sandwich for the linear block;
profile curvature for the sigma block){p_end}
{synopt:{cmd:e(sigma_se)}}profile-curvature SEs of sigma1..sigmaG{p_end}

{marker examples}{...}
{title:Examples}

{pstd}From raw shares, single nest, all regressors exogenous:{p_end}

{phang2}{cmd:. exactnl y x1 x2, nest(nestid) market(mkt) share(s_j) outside(s_0)}{p_end}

{pstd}Two nests, fixed effects, price instrumented by a cost shifter via
control function, a common crowding term, and clustered SEs:{p_end}

{phang2}{cmd:. exactnl y x1 x2 (price = costshifter), nest(nestid) market(mkt year) ///}{p_end}
{phang2}{cmd:      absorb(mkt year) share(s_j) outside(s_0) crowding("lnJ") cluster(mkt)}{p_end}

{pstd}Uncorrected estimates for comparison (Jacobian dropped):{p_end}

{phang2}{cmd:. exactnl y x1 x2 (price = costshifter), nest(nestid) market(mkt year) ///}{p_end}
{phang2}{cmd:      absorb(mkt year) share(s_j) outside(s_0) crowding("lnJ") cluster(mkt) noJACobian}{p_end}

{pstd}For data that carry a precomputed within-nest log share instead of raw
shares, replace {cmd:share()}+{cmd:outside()} with {cmd:wnshare()}.  See
{cmd:exactnl_example.do} for a runnable script, including the specification
used in the paper's Medicare Advantage application.{p_end}

{marker remarks}{...}
{title:Remarks}

{pstd}
Documentation, updates, and issue reporting are hosted at
{browse "https://github.com/lisa-george/exactnl"}.  Please file bug reports and
feature requests there.

{marker authors}{...}
{title:Authors}

{pstd}
Matthew J. Baker and Lisa M. George, Hunter College and the Graduate Center,
CUNY.

{pstd}
Website: {browse "https://github.com/lisa-george/exactnl"}

{marker references}{...}
{title:References}

{phang}
Ackerberg, D. A., and M. Rysman. 2005. Unobserved product differentiation in
discrete-choice models: Estimating price elasticities and welfare effects.
{it:RAND Journal of Economics} 36: 771-788.

{phang}
Baker, M. J., and L. M. George. 2026. Demand estimation with variable choice sets:
A likelihood correction for nested logit. SSRN Working Paper 7207938. Available at
{browse "https://ssrn.com/abstract=7207938"}.

{phang}
Berry, S. 1994. Estimating discrete-choice models of product differentiation.
{it:RAND Journal of Economics} 25: 242-262.

{phang}
Kodde, D. A., and F. C. Palm. 1986. Wald criteria for jointly testing equality
and inequality restrictions. {it:Econometrica} 54: 1243-1246.

{phang}
Self, S. G., and K.-Y. Liang. 1987. Asymptotic properties of maximum likelihood
estimators and likelihood ratio tests under nonstandard conditions.
{it:Journal of the American Statistical Association} 82: 605-610.
