*==============================================================================
* exactnl_cert_public.do  --  PUBLIC certification script for the exactnl command.
*
* This is the shippable, self-contained certification: it runs ONLY the simulated
* Monte Carlo check (CHECK 1 of the full private certification).  It uses NO
* Medicare Advantage data -- everything is generated from a fixed seed -- so it
* can be distributed and re-run by anyone who installs exactnl.
*
*   CHECK 1: Monte Carlo no-congestion DGP (sigma_true=0.5, G=1).  exactnl must
*            recover sigma and MATCH the paper's simulation-code profile estimate
*            on the same seed and the same data.
*
* The full three-check certification (which additionally reproduces the MA
* estimation-sample estimates and the noJACobian comparison) is kept private in
* exactnl_cert.do and requires the confidential estimation sample.
*
* This check uses raw shares and no price instrument, so it requires only
* exactnl itself plus reghdfe/ftools.
*==============================================================================
version 17
clear all
set more off

* ---- paths ---------------------------------------------------------------
* If exactnl is installed (net install exactnl), no path setup is needed and
* results are written to the current working directory. If running from a
* cloned copy of the repository, set EXACTNL_HOME to the package folder.
global EXACTNL_HOME "."
adopath ++ "$EXACTNL_HOME"

* results holder
tempname results
file open `results' using "$EXACTNL_HOME/exactnl_cert_public_results.txt", write replace text
file write `results' "exactnl PUBLIC certification results" _n
file write `results' "====================================" _n _n

*==============================================================================
* CHECK 1 -- Monte Carlo no-congestion DGP, sigma_true = 0.5
*   Single nest per market (G=1) so that exactnl's per-nest sigma is directly
*   and exactly comparable to the paper MC code's single common sigma.  Both use
*   the identical corrected objective (Gaussian LL + sum(J-1)ln(1-sigma)).
*==============================================================================
display _n(2) as result "==================== CHECK 1: Monte Carlo (sigma=0.5) ===================="

set seed 20260505
local M = 1000
local sigma_true = 0.5
local alpha_true = -1
local beta_true  = 1
local gam = 0

* MA-calibrated cumulative J probabilities (identical to 08_montecarlo_dgp.do)
local cp1  = 7622 / 53745
local cp2  = (7622 + 8325) / 53745
local cp3  = (7622 + 8325 + 7403) / 53745
local cp4  = (7622 + 8325 + 7403 + 6509) / 53745
local cp5  = (7622 + 8325 + 7403 + 6509 + 5105) / 53745
local cp6  = (7622 + 8325 + 7403 + 6509 + 5105 + 3714) / 53745
local cp7  = (7622 + 8325 + 7403 + 6509 + 5105 + 3714 + 2930) / 53745
local cp8  = (7622 + 8325 + 7403 + 6509 + 5105 + 3714 + 2930 + 2341) / 53745
local cp9  = (7622 + 8325 + 7403 + 6509 + 5105 + 3714 + 2930 + 2341 + 1927) / 53745
local cp10 = (7622 + 8325 + 7403 + 6509 + 5105 + 3714 + 2930 + 2341 + 1927 + 1538) / 53745
local cp11 = (7622 + 8325 + 7403 + 6509 + 5105 + 3714 + 2930 + 2341 + 1927 + 1538 + 1206) / 53745
local cp12 = (7622 + 8325 + 7403 + 6509 + 5105 + 3714 + 2930 + 2341 + 1927 + 1538 + 1206 + 977) / 53745
local cp13 = (7622 + 8325 + 7403 + 6509 + 5105 + 3714 + 2930 + 2341 + 1927 + 1538 + 1206 + 977 + 806) / 53745
local cp14 = (7622 + 8325 + 7403 + 6509 + 5105 + 3714 + 2930 + 2341 + 1927 + 1538 + 1206 + 977 + 806 + 676) / 53745

clear
set obs `M'
gen int market = _n
gen int nest = 1                       // single nest per market (G=1)
gen double u_draw = runiform()
gen int J_g = 20
replace J_g = 1  if u_draw <= `cp1'
replace J_g = 2  if u_draw > `cp1'  & u_draw <= `cp2'  & J_g == 20
replace J_g = 3  if u_draw > `cp2'  & u_draw <= `cp3'  & J_g == 20
replace J_g = 4  if u_draw > `cp3'  & u_draw <= `cp4'  & J_g == 20
replace J_g = 5  if u_draw > `cp4'  & u_draw <= `cp5'  & J_g == 20
replace J_g = 6  if u_draw > `cp5'  & u_draw <= `cp6'  & J_g == 20
replace J_g = 7  if u_draw > `cp6'  & u_draw <= `cp7'  & J_g == 20
replace J_g = 8  if u_draw > `cp7'  & u_draw <= `cp8'  & J_g == 20
replace J_g = 9  if u_draw > `cp8'  & u_draw <= `cp9'  & J_g == 20
replace J_g = 10 if u_draw > `cp9'  & u_draw <= `cp10' & J_g == 20
replace J_g = 11 if u_draw > `cp10' & u_draw <= `cp11' & J_g == 20
replace J_g = 12 if u_draw > `cp11' & u_draw <= `cp12' & J_g == 20
replace J_g = 13 if u_draw > `cp12' & u_draw <= `cp13' & J_g == 20
replace J_g = 14 if u_draw > `cp13' & u_draw <= `cp14' & J_g == 20
drop u_draw

expand J_g
bysort market nest: gen int product = _n
gen double x_j  = rnormal(0, 1)
gen double p_j  = rnormal(2, 1)
gen double xi_j = rnormal(0, 1)
gen double delta_j = `beta_true'*x_j + `alpha_true'*p_j + xi_j - `sigma_true'*`gam'*ln(J_g)
gen double exp_d_s = exp(delta_j / (1 - `sigma_true'))
bysort market nest: egen double D_g = total(exp_d_s)
gen double D_g_pow = D_g^(1 - `sigma_true')
bysort market nest: gen byte _fn = (_n == 1)
gen double _Dp_once = D_g_pow * _fn
bysort market: egen double denom_mkt = total(_Dp_once)
drop _fn _Dp_once
gen double s_jg = exp_d_s / D_g
gen double s_g  = D_g_pow / (1 + denom_mkt)
gen double s_j  = s_jg * s_g
gen double s_0  = 1 / (1 + denom_mkt)
gen double y        = ln(s_j) - ln(s_0)
gen double ln_within = ln(s_jg)

* --- paper simulation-code profile (single sigma; copy of 08_montecarlo_dgp.do)
bysort market nest: gen byte _fn2 = (_n == 1)
gen int _jm1 = (J_g - 1) * _fn2
summarize _jm1, meanonly
local sum_jm1 = r(sum)
drop _fn2 _jm1

local best_LL_c = -1e30
local best_sig_c = 0.5
forvalues sig_c = 5(5)90 {
    local sig = `sig_c' / 100
    gen double _z = y - `sig' * ln_within
    quietly reg _z x_j p_j
    local SSR = e(rss)
    local N_r = e(N)
    local LL_u = -`N_r'/2 * ln(`SSR'/`N_r') - `N_r'/2 - `N_r'/2 * ln(2*_pi)
    local LL_t = `LL_u' + `sum_jm1' * ln(1 - `sig')
    if `LL_t' > `best_LL_c' {
        local best_LL_c = `LL_t'
        local best_sig_c = `sig'
    }
    drop _z
}
local f_lo = max(0.02, `best_sig_c' - 0.05)
local f_hi = min(0.98, `best_sig_c' + 0.05)
local best_LL_f = `best_LL_c'
local best_sig_f = `best_sig_c'
local sig_f = `f_lo'
while `sig_f' <= `f_hi' + 0.0001 {
    gen double _z = y - `sig_f' * ln_within
    quietly reg _z x_j p_j
    local SSR = e(rss)
    local N_r = e(N)
    local LL_u = -`N_r'/2 * ln(`SSR'/`N_r') - `N_r'/2 - `N_r'/2 * ln(2*_pi)
    local LL_t = `LL_u' + `sum_jm1' * ln(1 - `sig_f')
    if `LL_t' > `best_LL_f' {
        local best_LL_f = `LL_t'
        local best_sig_f = `sig_f'
    }
    drop _z
    local sig_f = `sig_f' + 0.005
}
local sig_paper = `best_sig_f'

* --- exactnl on the SAME data (0.005 lattice to match the paper's grid) --------
exactnl y x_j p_j, nest(nest) market(market) share(s_j) outside(s_0) ///
    grid(0.005 0.995 0.005)
local sig_exnl = _b[sigma1]

* --- exactnl at the fine 0.001 grid (best point estimate) ----------------------
exactnl y x_j p_j, nest(nest) market(market) share(s_j) outside(s_0)
local sig_exnl_fine = _b[sigma1]

local d_match = abs(`sig_exnl' - `sig_paper')
local d_true  = abs(`sig_exnl_fine' - `sigma_true')

display as result "  sigma_true       = " %7.4f `sigma_true'
display as result "  paper MC profile = " %7.4f `sig_paper'
display as result "  exactnl (0.005)  = " %7.4f `sig_exnl'   "   |exactnl-paper| = " %8.5f `d_match'
display as result "  exactnl (0.001)  = " %7.4f `sig_exnl_fine' "   |exactnl-0.5|   = " %8.5f `d_true'

* agreement to one grid step (0.005) is the honest standard for grid estimators;
* d_match is typically 0 (same lattice point).
local pass1 = (`d_match' < 0.0051) & (`d_true' < 0.06)
file write `results' "CHECK 1 -- Monte Carlo no-congestion DGP (sigma_true=0.5, G=1, seed 20260505)" _n
file write `results' "  paper simulation-code profile sigma = " %7.4f (`sig_paper') _n
file write `results' "  exactnl sigma (same 0.005 grid)      = " %7.4f (`sig_exnl') "   |diff| = " %9.6f (`d_match') _n
file write `results' "  exactnl sigma (fine 0.001 grid)      = " %7.4f (`sig_exnl_fine') "   |sigma-0.5| = " %7.4f (`d_true') _n
if `pass1' file write `results' "  RESULT: PASS (exactnl matches paper code to grid tolerance; recovers 0.5)" _n _n
else       file write `results' "  RESULT: SEE REPORT (differences above)" _n _n

*------------------------------------------------------------------------------
file write `results' "SUMMARY (public): check1=" (`pass1') _n
file close `results'
type "$EXACTNL_HOME/exactnl_cert_public_results.txt"

display _n(2) as result "CHECK 1 pass = `pass1'"
display _n(2) "EXACTNL_CERT_PUBLIC_DONE_SENTINEL"
