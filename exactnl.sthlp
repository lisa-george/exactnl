{smcl}
{* *! version 1.0.0  09jul2026  Baker & George}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "nlogit" "help nlogit"}{...}
{viewerjumpto "Syntax" "exactnl##syntax"}{...}
{viewerjumpto "Description" "exactnl##description"}{...}
{viewerjumpto "Options" "exactnl##options"}{...}
{viewerjumpto "Stored results" "exactnl##results"}{...}
{viewerjumpto "Examples" "exactnl##examples"}{...}
{viewerjumpto "Remarks" "exactnl##remarks"}{...}
{viewerjumpto "Authors" "exactnl##authors"}{...}
{viewerjumpto "References" "exactnl##references"}{...}
{title:Title}

{phang}
{bf:exactnl} {hline 2} Exact-likelihood nested logit with the (J-1)ln(1-sigma)
Jacobian correction, fixed-effect absorption, and control-function price IV

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
{synopt:{opth market(varlist)}}market identifiers (e.g. {cmd:fipscode year}); a
nest is a {it:market x nest()} cell{p_end}

{syntab:Within-nest shares (supply ONE of the two)}
{synopt:{opth share(varname)}}plan-level market share s_j; combined with
{cmd:outside()} the command builds ln s_(j|g) and validates the shares{p_end}
{synopt:{opth out:side(varname)}}outside-good share s_0 (required with
{cmd:share()}){p_end}
{synopt:{opth wns:hare(varname)}}precomputed within-nest log share ln s_(j|g),
for data that do not carry raw shares (the MA estimation sample){p_end}

{syntab:Specification}
{synopt:{opth absorb(varlist)}}fixed effects to absorb, reghdfe-style (e.g.
{cmd:fipscode year}){p_end}
{synopt:{opt crow:ding(str)}}{cmd:lnJ} adds tau*ln(J_gm) to mean utility with a
common tau; {cmd:lnJ, bynest} adds a nest-specific tau_g{p_end}
{synopt:{opt noJAC:obian}}omit the Jacobian correction (for comparison; prints a
warning that sigma is biased upward){p_end}
{synopt:{opt sigma0(numlist)}}starting values for the golden-section refinement
when G>=3 (default 0.2 per nest){p_end}

{syntab:Inference}
{synopt:{opth cluster(varname)}}cluster-robust (sandwich) SEs for the linear
parameters{p_end}
{synopt:{opt boot:strap(reps)}}county block-bootstrap SEs for sigma_g, tau and
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

{marker description}{...}
{title:Description}

{pstd}
{cmd:exactnl} estimates a two-level nested logit demand model by the exact
(corrected) profile-likelihood method rather than by the Berry (1994) linear IV
regression of the log share-ratio on the within-nest log share.  Berry's
instrumental-variables estimator treats ln s_(j|g) as a single endogenous
regressor and instruments it with the count of rival plans in the nest.  When
nest size is correlated with the unobserved mean utility (congestion, crowding),
that count instrument is invalid and the IV estimate of the nesting parameter is
biased.  The likelihood estimator writes the log likelihood of the observed
shares directly and includes the change-of-variables Jacobian of the nested-logit
share map.  For nest g in market m containing J_gm plans, that Jacobian
contributes {bf:(J_gm - 1) ln(1 - sigma_g)} to the log likelihood, which is the
term the linear IV regression discards.  Concentrating the linear parameters out
by fixed-effect partialling, the profile log likelihood is

{p 8 8 2}
LL(sigma) = -N/2 ln(SSR(sigma)/N) + sum_g [ sum_m (J_gm - 1) ] ln(1 - sigma_g),

{pstd}
where SSR(sigma) is a quadratic form in sigma built once from the FE-partialled
residuals of {it:depvar} and of the within-nest log shares interacted with each
nest.  {cmd:exactnl} maximises LL over sigma in [lo,hi]^G (an exhaustive grid for
G<=2, matching the paper's implementation exactly; coordinate golden-section
ascent for G>=3), then recovers the linear parameters (alpha on price, beta on
{it:indepvars}, and tau on the crowding term) by one final FE regression at the
optimum.  This reproduces the corrected-likelihood column of Table 4 of the MA
nesting paper.

{pstd}
Price endogeneity is handled by a control function: the first stage regresses
{it:endogvar} on the excluded instrument(s), the exogenous regressors, the
crowding term and the absorbed fixed effects (via {help reghdfe}), and the
first-stage residual is added to the structural equation.  {cmd:reghdfe} is
invoked directly for the first stage (not {cmd:ivreghdfe}) so that the residual
is available as a generated regressor.  Standard errors on the linear parameters
are the {cmd:cluster()} sandwich from the final regression; the sigma_g standard
errors are profile-curvature (observed-information) approximations, and a cluster
block bootstrap ({cmd:bootstrap()}) is recommended for inference on sigma_g.  The
command also reports the boundary likelihood-ratio test of H0: all sigma_g = 0
against the one-sided alternative.  The printed 5% critical value is the
Kodde-Palm (1986) least-favourable upper bound (2.71 for G=1, 5.14 for G=2),
which is valid regardless of the dependence among the sigma_g estimates; the
exact chi-bar-squared mixture quantile (valid under independence of the sigma_g
estimates) is stored in {cmd:e(LR_crit5_exact)} (see
{help exactnl##results:Stored results}).

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
{opth share(varname)} / {opth outside(varname)} supply the raw plan share s_j and
outside share s_0.  The command builds ln s_(j|g) = ln(s_j) - ln(sum of s_j in the
nest) and validates that shares lie in (0,1) and that within-market share totals
do not exceed 1.

{phang}
{opth wnshare(varname)} supplies a precomputed within-nest log share ln s_(j|g)
directly.  Use this for data (such as the MA estimation sample) that store the
within-nest log share but not the raw shares.  Exactly one of {cmd:wnshare()} or
{cmd:share()}+{cmd:outside()} must be given.

{phang}
{opth absorb(varlist)} lists the fixed effects to absorb ({help reghdfe} syntax).

{phang}
{opt crowding(str)} adds a nest-size (crowding) term to mean utility.
{cmd:crowding(lnJ)} adds tau*ln(J_gm) with a single common tau;
{cmd:crowding(lnJ, bynest)} adds a separate tau_g for each nest.

{phang}
{opt noJACobian} drops the (J_gm-1)ln(1-sigma_g) term, reducing the objective to
the Gaussian sum of squares.  The resulting sigma_g are the uncorrected
(upward-biased) estimates; a warning is printed.

{phang}
{opt sigma0(numlist)} sets starting values for the G>=3 golden-section refinement.

{phang}
{opth cluster(varname)} requests cluster-robust SEs for the linear parameters.

{phang}
{opt bootstrap(reps)} runs a county (cluster) block bootstrap of the full
estimator, resampling clusters, rebuilding the nest counts on the bootstrap
replica id, and refitting; it reports bootstrap SEs for sigma_g, tau and the
linear parameters.  Requires {cmd:cluster()} to be the leading {cmd:market()}
variable.

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

{pstd}MA application, preferred spec (two nests HMO/PPO, common crowding tau,
county+year FE, county-clustered SEs, price instrumented by the CMS benchmark):{p_end}

{phang2}{cmd:. gen byte gap_coverage = (extragap=="Y")}{p_end}
{phang2}{cmd:. keep if share_flag==0}{p_end}
{phang2}{cmd:. exactnl depvar_logit ddrugdeduct gap_coverage (premium_partc = benchmark_no_bonus), ///}{p_end}
{phang2}{cmd:      nest(plan_group) market(fipscode year) absorb(fipscode year) ///}{p_end}
{phang2}{cmd:      wnshare(ln_within_nest_share) crowding("lnJ") cluster(fipscode)}{p_end}

{pstd}Uncorrected (Berry-style) sigma for comparison:{p_end}

{phang2}{cmd:. exactnl depvar_logit ddrugdeduct gap_coverage (premium_partc = benchmark_no_bonus), ///}{p_end}
{phang2}{cmd:      nest(plan_group) market(fipscode year) absorb(fipscode year) ///}{p_end}
{phang2}{cmd:      wnshare(ln_within_nest_share) crowding("lnJ") cluster(fipscode) noJACobian}{p_end}

{pstd}From raw shares, single nest, no price instrument:{p_end}

{phang2}{cmd:. exactnl y x1 x2, nest(one) market(mkt) share(s_j) outside(s_0)}{p_end}

{pstd}See {cmd:exactnl_example.do} for a runnable script.{p_end}

{marker remarks}{...}
{title:Remarks}

{pstd}
Documentation, updates, and issue reporting are hosted at
{browse "https://github.com/lisa-george/exactnl"}.  Please file bug reports and
feature requests there.

{marker authors}{...}
{title:Authors}

{pstd}
Packaged from the corrected-likelihood estimation core of the Medicare Advantage
nesting project (Baker and George).  Correspondence: lisa.m.george@gmail.com.

{pstd}
Website: {browse "https://github.com/lisa-george/exactnl"}

{marker references}{...}
{title:References}

{phang}
Baker and George. 2024. Nesting, crowding, and the identification of the
nested-logit correlation parameter in Medicare Advantage. Working paper.

{phang}
Berry, S. 1994. Estimating discrete-choice models of product differentiation.
{it:RAND Journal of Economics} 25: 242-262.

{phang}
Kodde, D. A., and F. C. Palm. 1986. Wald criteria for jointly testing equality
and inequality restrictions. {it:Econometrica} 54: 1243-1248.

{phang}
Self, S. G., and K.-Y. Liang. 1987. Asymptotic properties of maximum likelihood
estimators and likelihood ratio tests under nonstandard conditions.
{it:Journal of the American Statistical Association} 82: 605-610.
