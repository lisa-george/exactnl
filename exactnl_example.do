*==============================================================================
* exactnl_example.do
* Worked examples for the exactnl command (exact maximum likelihood estimation
* of the share-form nested logit).
*
* Section A is self-contained: it simulates nested logit data and runs the
* command. It requires exactnl (net install exactnl) plus reghdfe and ftools
* (ssc install reghdfe; ssc install ftools). Section B shows the specification
* used in the paper's Medicare Advantage application; it requires the
* confidential estimation sample and is provided for reference.
*==============================================================================
version 17
clear all
set more off
set seed 20260710

*==============================================================================
* SECTION A -- simulated data (runs anywhere)
*==============================================================================

*--- generate a two-nest nested logit dataset with known parameters ----------
* True values: sigma = 0.5 (both nests), coefficient on x = 1, coefficient on
* price = -1. Price is endogenous (correlated with xi); z is a cost-shifter
* instrument.
local M = 500                          // markets
set obs `M'
gen long mkt = _n
gen int  J   = 2 + floor(9*runiform())  // 2-10 products per market
expand J
bysort mkt: gen int prod = _n
gen byte nestid = 1 + mod(prod, 2)      // alternate products across two nests
gen double xi = 0.3*rnormal()
gen double x  = rnormal()
gen double z  = rnormal()
gen double p  = 1 + 0.5*z + 0.5*xi + 0.2*rnormal()   // endogenous price
gen double delta = 1 + x - p + xi

* nested logit shares from the closed form, sigma_true = 0.5
local sig = 0.5
gen double ed = exp(delta/(1-`sig'))
bysort mkt nestid: egen double Dg = total(ed)
gen double snum = ed * Dg^(-`sig')                    // = ed * Dg^{1-sig} / Dg
bysort mkt nestid: gen byte firstnest = (_n==1)
egen double sumD = total(cond(firstnest, Dg^(1-`sig'), 0)), by(mkt)
gen double s_j = snum / (1 + sumD)
gen double s_0 = 1 / (1 + sumD)
gen double y   = ln(s_j) - ln(s_0)

*--- Example 1: corrected likelihood, price by control function --------------
* Expect: sigma1 and sigma2 near 0.5, the coefficient on p near -1, and the
* coefficient on x near 1.
exactnl y x (p = z), nest(nestid) market(mkt) share(s_j) outside(s_0) ///
    cluster(mkt)
ereturn list

*--- Example 2: uncorrected comparison (Jacobian dropped) ---------------------
* Expect sigma biased upward; the command prints a warning.
exactnl y x (p = z), nest(nestid) market(mkt) share(s_j) outside(s_0) ///
    cluster(mkt) noJACobian

*--- Example 3: add a common crowding term tau*ln(J_gm) -----------------------
* The DGP has no crowding, so expect tau near 0 and sigma unchanged.
exactnl y x (p = z), nest(nestid) market(mkt) share(s_j) outside(s_0) ///
    crowding("lnJ") cluster(mkt)

*==============================================================================
* SECTION B -- full real-data syntax: the specification from the paper's
* Medicare Advantage application (Baker and George, working paper).
* Requires the paper's estimation sample; shown for syntax reference.
* bootstrap(200) reproduces the county block-bootstrap standard errors the
* paper reports for sigma_g and tau. It refits the full model 200 times and
* is compute-intensive on the full sample; drop it for a quick point estimate.
*==============================================================================
/*
exactnl depvar_logit ddrugdeduct gap_coverage ///
        (premium_partc = benchmark_no_bonus), ///
    nest(plan_group) market(fipscode year) absorb(fipscode year) ///
    wnshare(ln_within_nest_share) crowding("lnJ") ///
    cluster(fipscode) bootstrap(200)
*/
