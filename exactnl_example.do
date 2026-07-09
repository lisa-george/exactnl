*==============================================================================
* exactnl_example.do
* Worked examples for the exactnl command (exact-likelihood nested logit with
* the (J-1)ln(1-sigma) Jacobian correction).
*
* Set the two globals below, then run.  The command file exactnl.ado must be on
* the adopath (the adopath line does this if you keep this .do next to the .ado).
*==============================================================================
version 17
clear all
set more off

* --- point Stata at the exactnl package (folder containing exactnl.ado) --------
adopath ++ "D:/Dropbox/NML/nml_lisa_matt_shared/exactnl"

* --- the frozen MA estimation sample -------------------------------------------
global SAMPLE "D:/Dropbox/NML/nml_lisa_matt_shared/data_medicare/clean/ma_nml_estimation_sample.dta"

use "$SAMPLE", clear
gen byte gap_coverage = (extragap=="Y" | extragap=="Yes")
keep if share_flag==0        // drop problematic-share markets (paper filter)

*------------------------------------------------------------------------------
* Example 1.  Preferred specification, Table 4 corrected-likelihood column.
*   two nests (HMO=1, PPO=2), common crowding term tau*ln(J), county+year FE,
*   price (premium_partc) instrumented by the CMS benchmark via control function,
*   county-clustered SEs.  Reproduces sigma_HMO=0.136, sigma_PPO=0.117,
*   tau=-0.052, alpha=-0.0336.
*------------------------------------------------------------------------------
exactnl depvar_logit ddrugdeduct gap_coverage i.year ///
        (premium_partc = benchmark_no_bonus), ///
    nest(plan_group) market(fipscode year) absorb(fipscode) ///
    wnshare(ln_within_nest_share) crowding("lnJ") cluster(fipscode)

* stored results
ereturn list
matrix list e(b)

*------------------------------------------------------------------------------
* Example 2.  Uncorrected (Berry-style) sigma for comparison: same spec, but
*   drop the Jacobian.  sigma is biased upward.
*------------------------------------------------------------------------------
exactnl depvar_logit ddrugdeduct gap_coverage i.year ///
        (premium_partc = benchmark_no_bonus), ///
    nest(plan_group) market(fipscode year) absorb(fipscode) ///
    wnshare(ln_within_nest_share) crowding("lnJ") cluster(fipscode) noJACobian

*------------------------------------------------------------------------------
* Example 3.  Nest-specific crowding (separate tau_HMO, tau_PPO).
*------------------------------------------------------------------------------
exactnl depvar_logit ddrugdeduct gap_coverage i.year ///
        (premium_partc = benchmark_no_bonus), ///
    nest(plan_group) market(fipscode year) absorb(fipscode) ///
    wnshare(ln_within_nest_share) crowding("lnJ, bynest") cluster(fipscode)

*------------------------------------------------------------------------------
* Example 4.  No crowding term (plain corrected nested logit).
*------------------------------------------------------------------------------
exactnl depvar_logit ddrugdeduct gap_coverage i.year ///
        (premium_partc = benchmark_no_bonus), ///
    nest(plan_group) market(fipscode year) absorb(fipscode) ///
    wnshare(ln_within_nest_share) cluster(fipscode)

*------------------------------------------------------------------------------
* Example 5.  Block-bootstrap SEs for sigma_g and tau (county blocks).
*   (Commented out: this reruns the estimator `reps' times.)
*------------------------------------------------------------------------------
* exactnl depvar_logit ddrugdeduct gap_coverage i.year ///
*         (premium_partc = benchmark_no_bonus), ///
*     nest(plan_group) market(fipscode year) absorb(fipscode) ///
*     wnshare(ln_within_nest_share) crowding("lnJ") cluster(fipscode) ///
*     bootstrap(200)

display _n(2) "exactnl_example.do complete."
