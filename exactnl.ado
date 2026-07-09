*! version 1.0.0  09jul2026  Baker & George
*! exactnl : exact-likelihood nested logit with the (J-1)ln(1-sigma) Jacobian
*!           correction, FE absorption, and control-function price IV.
*! Estimation core extracted from the corrected-likelihood two-sigma do-file
*! behind Table 4 of Baker & George (MA nesting paper); see exactnl.sthlp.
*
* Model (per market m, nest g, plan j):
*   ln(s_j) - ln(s_0) = X'beta - alpha*price + sigma_g * ln(s_{j|g,m})
*                       [ + tau * ln(J_gm) | + tau_g * ln(J_gm) ] + xi
* Price is handled by a control function: v_hat = residual of price on the
* excluded instrument(s) + exogenous regressors + FE (reghdfe first stage),
* added to the structural equation.  Given (sigma_1..sigma_G) the equation is
* linear, so the profile (concentrated) log likelihood is
*   LL(sigma) = -N/2 ln(SSR(sigma)/N) + sum_g (J_gm-1 summed over m) ln(1-sigma_g)
* maximised over sigma in [lo,hi]^G.  SSR(sigma) is a quadratic form in sigma
* built once from FE-partialled residuals of the depvar and the within-nest
* log shares (algebraically identical to re-running reghdfe at every grid node).

program define exactnl, eclass
    version 17

    * ---- separate command body (may contain "(endog = inst)") from options ----
    gettoken body opts : 0, parse(",")
    if `"`opts'"' != "" {
        * strip the leading comma
        local opts = substr(`"`opts'"', 2, .)
    }

    * ---- detect and remove the negatable noJACobian flag manually.  It takes no
    * argument, so stripping it before syntax() is safe and avoids parser quirks.
    local nojac = 0
    if ustrregexm(`"`opts'"', "(?i)\bnojac") {
        local nojac = 1
        local opts = ustrregexra(`"`opts'"', "(?i)\bnojac\w*", " ")
    }

    * ---- pull the control-function IV group "(endog = instruments)" ----
    local endog ""
    local inst  ""
    if regexm(`"`body'"', "\(([^)]*)\)") {
        local ivgrp = regexs(1)
        local body  = regexr(`"`body'"', "\([^)]*\)", " ")
        gettoken endog inst : ivgrp, parse("=")
        local inst = subinstr(`"`inst'"', "=", "", 1)
        local endog = strtrim(`"`endog'"')
        local inst  = strtrim(`"`inst'"')
    }

    * ---- depvar + exogenous indepvars from what remains of the body ----
    local body = strtrim(stritrim(`"`body'"'))
    gettoken depvar indepvars : body
    local indepvars = strtrim(`"`indepvars'"')

    * ---- options ----
    local 0 `", `opts'"'
    syntax , NEST(varname) MARKET(varlist) ///
        [ SHARE(varname) OUTside(varname) WNShare(varname) ///
          ABSORB(varlist) CROWding(string) SIGMA0(numlist) ///
          CLUSTER(varname) GRID(numlist) ///
          BOOTstrap(integer 0) TOLerance(real 1e-8) ]

    * validate depvar
    confirm numeric variable `depvar'

    * ---- within-nest log share: supplied directly, or built from raw shares --
    if "`wnshare'" != "" & "`share'" != "" {
        di as error "specify either wnshare() OR share()+outside(), not both"
        exit 198
    }
    if "`wnshare'" == "" & "`share'" == "" {
        di as error "you must supply the within-nest log share via wnshare(), "
        di as error "or raw shares via share() and outside()"
        exit 198
    }

    * ---- crowding option ----
    local crowd_on  = 0
    local crowd_by  = 0
    if `"`crowding'"' != "" {
        local ctok = lower(subinstr(strtrim(`"`crowding'"'), " ", "", .))
        if strpos("`ctok'", "lnj") == 0 {
            di as error "crowding() must be lnJ  or  lnJ, bynest"
            exit 198
        }
        local crowd_on = 1
        if strpos("`ctok'", "bynest") > 0 local crowd_by = 1
    }

    * ---- grid ----
    if "`grid'" == "" {
        local glo = 0.001
        local ghi = 0.999
        local gst = 0.001
    }
    else {
        tokenize "`grid'"
        local glo `1'
        local ghi `2'
        local gst `3'
        if "`gst'" == "" {
            di as error "grid() takes three numbers: lo hi step"
            exit 198
        }
    }

    * Structural sample: everything the mean-utility equation needs EXCEPT the
    * excluded instrument.  Nest counts J_gm and the Jacobian weights are market
    * primitives (how many plans compete in each nest) and must be built on this
    * sample -- a plan with a missing instrument is still a competitor and is
    * still counted in J_gm.  This matches the reference esample, which never
    * conditions on the instrument.  The instrument sample is imposed afterwards,
    * only for the first stage and the Gram.
    marksample touse, novarlist
    markout `touse' `depvar' `indepvars' `endog' `nest' `market' ///
        `wnshare' `share' `outside' `absorb' `cluster'

    * ============================ build regressors ==========================
    tempvar wns lnJ Jcell tagcell

    * within-nest log share
    if "`wnshare'" != "" {
        gen double `wns' = `wnshare' if `touse'
    }
    else {
        * validate raw shares
        quietly count if `touse' & (`share'<=0 | `share'>=1)
        if r(N) > 0 {
            di as error "share() has `r(N)' obs outside (0,1)"
            exit 459
        }
        quietly count if `touse' & (`outside'<=0 | `outside'>=1)
        if r(N) > 0 {
            di as error "outside() has `r(N)' obs outside (0,1)"
            exit 459
        }
        tempvar mktsum
        bysort `touse' `market' (`nest'): egen double `mktsum' = total(`share') if `touse'
        quietly count if `touse' & `mktsum' > 1 + 1e-6 & !missing(`mktsum')
        if r(N) > 0 {
            di as error "within-market share totals exceed 1 in `r(N)' markets"
            exit 459
        }
        tempvar nestsum
        bysort `touse' `market' `nest': egen double `nestsum' = total(`share') if `touse'
        gen double `wns' = ln(`share') - ln(`nestsum') if `touse'
    }

    * nest counts J_gm and crowding term
    bysort `touse' `market' `nest': gen long `Jcell' = _N if `touse'
    gen double `lnJ' = ln(`Jcell') if `touse'
    egen byte `tagcell' = tag(`market' `nest') if `touse'

    * nest levels -> per-nest within-share and crowding variables + Jacobian sums
    quietly levelsof `nest' if `touse', local(nlevels)
    local G : word count `nlevels'
    if `G' < 1 {
        di as error "nest() has no non-missing levels on the estimation sample"
        exit 459
    }

    local lwvars ""
    local crowdvars ""
    local k = 0
    tempname SUMv
    matrix `SUMv' = J(1, `G', 0)
    foreach lv of local nlevels {
        local ++k
        tempvar lw`k'
        gen double `lw`k'' = `wns' * (`nest' == `lv') if `touse'
        local lwvars `lwvars' `lw`k''
        * Jacobian sum for this nest: sum over cells of (J_gm - 1) = sum(J)-Ncells
        quietly total `Jcell' if `touse' & `tagcell' & `nest'==`lv'
        tempname tb
        matrix `tb' = e(b)
        quietly count if `touse' & `tagcell' & `nest'==`lv'
        local ncell_g = r(N)
        matrix `SUMv'[1,`k'] = `tb'[1,1] - `ncell_g'
        if `crowd_on' & `crowd_by' {
            tempvar cj`k'
            gen double `cj`k'' = `lnJ' * (`nest'==`lv') if `touse'
            local crowdvars `crowdvars' `cj`k''
        }
    }
    if `crowd_on' & !`crowd_by' {
        local crowdvars `lnJ'
    }

    * Now impose the instrument sample for the first stage and the Gram (J_gm and
    * SUMv above were built on the wider structural sample, as in the reference).
    if "`inst'" != "" {
        markout `touse' `inst'
    }

    * ============================ first stage (CF) ==========================
    tempvar vhat
    if "`endog'" != "" {
        * control-function first stage.  reghdfe is used explicitly (not
        * ivreghdfe) so the second-stage residual enters as a generated regressor.
        if "`absorb'" != "" {
            quietly reghdfe `endog' `inst' `indepvars' `crowdvars' if `touse', ///
                absorb(`absorb') residuals(`vhat')
        }
        else {
            quietly reghdfe `endog' `inst' `indepvars' `crowdvars' if `touse', ///
                noabsorb residuals(`vhat')
        }
        local rhs `indepvars' `endog' `crowdvars' `vhat'
    }
    else {
        * no endogenous regressor: exogenous OLS path
        local rhs `indepvars' `crowdvars'
    }

    * ============ FE-partialled residuals of depvar and each lw_g ============
    tempvar e0
    _exactnl_resid `depvar' `touse' "`rhs'" "`absorb'" `e0'
    local reslist `e0'
    local k = 0
    foreach lw of local lwvars {
        local ++k
        tempvar e`k'
        _exactnl_resid `lw' `touse' "`rhs'" "`absorb'" `e`k''
        local reslist `reslist' `e`k''
    }

    * common sample across all residualisations
    markout `touse' `reslist'
    quietly count if `touse'
    local N = r(N)
    quietly egen long `Jcell'2 = group(`market') if `touse'
    quietly summarize `Jcell'2, meanonly
    local Nmk = r(max)
    drop `Jcell'2

    * ============================ optimise sigma ============================
    local jac = (`nojac' == 0)
    tempname LLbest LL0 LR
    capture matrix drop __exactnl_sig __exactnl_sigse
    scalar __exactnl_LLbest = .
    scalar __exactnl_LL0    = .
    mata: exactnl_fit("`reslist'", "`touse'", "`SUMv'", `N', `glo', `ghi', ///
        `gst', `jac', `G')

    * sigma point estimates + curvature SEs returned in the __exactnl_* objects
    tempname sig sigse
    matrix `sig'   = __exactnl_sig
    matrix `sigse' = __exactnl_sigse
    scalar `LLbest' = __exactnl_LLbest
    scalar `LL0'    = __exactnl_LL0
    scalar `LR'     = 2*(`LLbest' - `LL0')
    scalar drop __exactnl_LLbest __exactnl_LL0
    matrix drop __exactnl_sig __exactnl_sigse

    * ======================= structural linear params =======================
    tempvar zstar
    quietly gen double `zstar' = `depvar' if `touse'
    local k = 0
    foreach lw of local lwvars {
        local ++k
        quietly replace `zstar' = `zstar' - `sig'[1,`k']*`lw' if `touse'
    }
    if "`cluster'" != "" local clopt cluster(`cluster')
    if "`absorb'" != "" {
        quietly reghdfe `zstar' `rhs' if `touse', absorb(`absorb') `clopt'
    }
    else {
        quietly reghdfe `zstar' `rhs' if `touse', noabsorb `clopt'
    }
    tempname bl Vl
    matrix `bl' = e(b)
    matrix `Vl' = e(V)

    * structural linear params to report (by NAME, so i.year factor terms and
    * the control-function vhat / _cons are excluded cleanly).  keepvars = the
    * actual (temp)var to pull from e(b); keeplab = its display name.
    local keepvars ""
    local keeplab  ""
    foreach v of local indepvars {
        if strpos("`v'",".")==0 & strpos("`v'","#")==0 {
            local keepvars `keepvars' `v'
            local keeplab  `keeplab'  `v'
        }
    }
    if "`endog'" != "" {
        local keepvars `keepvars' `endog'
        local keeplab  `keeplab'  `endog'
    }
    if `crowd_on' & !`crowd_by' {
        local keepvars `keepvars' `lnJ'
        local keeplab  `keeplab'  tau
    }
    else if `crowd_on' & `crowd_by' {
        local kk = 0
        foreach cv of local crowdvars {
            local ++kk
            local keepvars `keepvars' `cv'
            local keeplab  `keeplab'  tau`kk'
        }
    }

    * ============================ assemble e(b),e(V) ========================
    * order: linear params, then sigma_1..sigma_G
    local p : word count `keeplab'
    local q = `G'
    tempname b V
    matrix `b' = J(1, `p'+`q', 0)
    matrix `V' = J(`p'+`q', `p'+`q', 0)
    local cn ""
    forvalues i = 1/`p' {
        local vi : word `i' of `keepvars'
        local ni : word `i' of `keeplab'
        local ci = colnumb(`bl', "`vi'")
        matrix `b'[1,`i'] = `bl'[1,`ci']
        local cn `cn' `ni'
        forvalues jj = 1/`p' {
            local vj : word `jj' of `keepvars'
            local cj = colnumb(`bl', "`vj'")
            matrix `V'[`i',`jj'] = `Vl'[`ci',`cj']
        }
    }
    forvalues g = 1/`q' {
        local col = `p' + `g'
        matrix `b'[1,`col'] = `sig'[1,`g']
        matrix `V'[`col',`col'] = `sigse'[1,`g']^2
        local cn `cn' sigma`g'
    }
    matrix colnames `b' = `cn'
    matrix colnames `V' = `cn'
    matrix rownames `V' = `cn'

    ereturn post `b' `V', esample(`touse') depname(`depvar')
    ereturn scalar N        = `N'
    ereturn scalar N_markets = `Nmk'
    ereturn scalar G        = `G'
    ereturn scalar LL       = `LLbest'
    ereturn scalar LL_sigma0 = `LL0'
    ereturn scalar LR       = `LR'

    * ---- boundary LR 5% critical values (H0: all sigma_g = 0, one-sided) ----
    * (a) Kodde-Palm (1986) least-favourable upper bound: c solving
    *     0.5*chi2tail(G-1,c) + 0.5*chi2tail(G,c) = 0.05  (2.71 for G=1;
    *     5.14 for G=2).  Conservative, valid regardless of the dependence
    *     among the sigma_g estimates -- this is the value printed.
    local kplo = 0
    local kphi = 60
    forvalues it = 1/200 {
        local cmid = (`kplo'+`kphi')/2
        local t1 = 0
        if `G'-1 >= 1 local t1 = chi2tail(`G'-1,`cmid')
        local sv = 0.5*`t1' + 0.5*chi2tail(`G',`cmid')
        if (`sv' > 0.05) local kplo = `cmid'
        else             local kphi = `cmid'
    }
    ereturn scalar LR_crit5 = (`kplo'+`kphi')/2

    * (b) exact chi-bar-squared mixture quantile  sum_k C(G,k)/2^G chi2_k
    *     (2.71 for G=1; ~4.23 for G=2).  Exact only when the sigma_g
    *     estimates are independent; stored for reference as LR_crit5_exact.
    local clo0 = 0
    local chi0 = 60
    forvalues it = 1/200 {
        local cmid = (`clo0'+`chi0')/2
        local sv = 0
        forvalues kk = 1/`G' {
            local wk = comb(`G',`kk')/(2^`G')
            local sv = `sv' + `wk'*chi2tail(`kk',`cmid')
        }
        if (`sv' > 0.05) local clo0 = `cmid'
        else             local chi0 = `cmid'
    }
    ereturn scalar LR_crit5_exact = (`clo0'+`chi0')/2
    ereturn local  crowding `"`crowding'"'
    ereturn local  cmd  "exactnl"
    ereturn local  depvar "`depvar'"

    * store sigma SEs as a separate vector for the display table
    tempname sesig
    matrix `sesig' = `sigse'
    ereturn matrix sigma_se = `sesig'

    * ============================== display =================================
    _exactnl_display, jac(`jac') g(`G') linnames(`linnames') ///
        crowd(`crowd_on') crowdby(`crowd_by')

    if `bootstrap' > 0 {
        di as text _n "{p}Running county block bootstrap of sigma_g, tau, and " ///
            "the linear parameters (`bootstrap' reps)...{p_end}"
        _exactnl_boot, reps(`bootstrap') cluster(`cluster') ///
            depvar(`depvar') indepvars(`indepvars') endog(`endog') inst(`inst') ///
            nest(`nest') market(`market') absorb(`absorb') wns(`wns') ///
            touse(`touse') glo(`glo') ghi(`ghi') gst(`gst') jac(`jac') ///
            crowdon(`crowd_on') crowdby(`crowd_by')
    }
end

* ---------------------------------------------------------------------------
* helper: FE-partial one variable on rhs, return residual in `res'
program define _exactnl_resid
    args y touse rhs absorb res
    if "`absorb'" != "" {
        quietly reghdfe `y' `rhs' if `touse', absorb(`absorb') residuals(`res')
    }
    else {
        quietly reghdfe `y' `rhs' if `touse', noabsorb residuals(`res')
    }
end

* ---------------------------------------------------------------------------
* helper: format the results table
program define _exactnl_display
    syntax , jac(integer) g(integer) [ linnames(string) crowd(integer 0) ///
        crowdby(integer 0) ]
    di as text _n "{hline 78}"
    if `jac' {
        di as text "Exact-likelihood nested logit (Jacobian-corrected)"
    }
    else {
        di as text "Nested logit, UNCORRECTED likelihood (noJACobian)"
    }
    di as text "{hline 78}"
    di as text %-24s "Parameter" _col(30) %12s "Coef." _col(44) %12s "Std. Err." ///
        _col(58) "  z"
    di as text "{hline 78}"
    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    local names : colnames `b'
    local i = 0
    foreach nm of local names {
        local ++i
        local coef = `b'[1,`i']
        local se   = sqrt(`V'[`i',`i'])
        local z    = `coef'/`se'
        if substr("`nm'",1,5) == "sigma" {
            di as result %-24s "`nm'" _col(30) %12.4f `coef' _col(44) ///
                %12.4f `se' _col(58) %8.2f `z' as text "  (profile)"
        }
        else {
            di as result %-24s "`nm'" _col(30) %12.4f `coef' _col(44) ///
                %12.4f `se' _col(58) %8.2f `z'
        }
    }
    di as text "{hline 78}"
    di as text "N = " as result e(N) as text "   markets = " as result ///
        e(N_markets) as text "   nests G = " as result e(G)
    if `jac' {
        di as text "Boundary LR test  H0: all sigma_g = 0"
        di as text "   LR = " as result %8.2f e(LR) as text ///
            "   Kodde-Palm 5% bound = " as result %6.2f e(LR_crit5) ///
            as text "   ==> " _c
        if e(LR) > e(LR_crit5) di as result "reject H0 (nesting present)"
        else                   di as result "fail to reject H0"
    }
    di as text "{p 0 3 2}Note: sigma_g std. errors are profile-curvature " ///
        "(observed-information) approximations.  Cluster block-bootstrap SEs " ///
        "are recommended for inference on sigma_g; use option " ///
        "{bf:bootstrap(reps) cluster(varname)}.{p_end}"
    if `jac' == 0 {
        di as text "{p 0 3 2}WARNING: noJACobian omits the (J-1)ln(1-sigma) " ///
            "correction; sigma_g are biased upward and are for comparison " ///
            "only.{p_end}"
    }
end

* ---------------------------------------------------------------------------
* helper: county block bootstrap of sigma_g, tau, and linear params.
* Canonical design (matches the block bootstrap in the tab4 do-file): resample
* clusters, rebuild the nest counts J_gm on the bootstrap replica id, refit.
* Assumes cluster() is the leading variable of market() (the county in the MA
* application, market(fipscode year) cluster(fipscode)).
program define _exactnl_boot, rclass
    syntax , reps(integer) [ cluster(string) depvar(string) indepvars(string) ///
        endog(string) inst(string) nest(string) market(string) absorb(string) ///
        wns(string) touse(string) glo(real 0.001) ghi(real 0.999) ///
        gst(real 0.001) jac(integer 1) crowdon(integer 0) crowdby(integer 0) ]

    if "`cluster'" == "" {
        di as error "bootstrap() requires cluster()"
        exit 198
    }
    local othermkt : list market - cluster
    local crwopt ""
    if `crowdon' & `crowdby'  local crwopt crowding("lnJ, bynest")
    else if `crowdon'         local crwopt crowding("lnJ")
    local jacopt ""
    if `jac' == 0 local jacopt noJACobian
    local ivopt ""
    if "`endog'" != "" local ivopt (`endog' = `inst')

    preserve
    quietly keep if `touse'
    tempfile snap
    quietly save `snap'

    local B = `reps'
    local first = 1
    local names ""
    forvalues b = 1/`B' {
        quietly use `snap', clear
        bsample, cluster(`cluster') idcluster(__bid)
        capture noisily exactnl `depvar' `indepvars' `ivopt', ///
            nest(`nest') market(__bid `othermkt') absorb(`absorb') ///
            wnshare(`wns') `crwopt' `jacopt' grid(`glo' `ghi' `gst')
        if _rc == 0 {
            matrix __brow = e(b)
            if `first' {
                matrix __boot_acc = __brow
                local names : colnames __brow
                local first = 0
            }
            else {
                matrix __boot_acc = __boot_acc \ __brow
            }
        }
    }
    restore

    local ncol : word count `names'
    di as text _n "{hline 60}"
    di as text "Block-bootstrap SEs  (reps completed: " ///
        as result rowsof(__boot_acc) as text ", cluster " as result "`cluster'" as text ")"
    di as text "{hline 60}"
    di as text %-24s "Parameter" _col(30) %14s "Boot SE"
    preserve
    forvalues j = 1/`ncol' {
        local nm : word `j' of `names'
        tempname col
        matrix `col' = __boot_acc[1...,`j']
        clear
        quietly svmat double `col'
        quietly summarize `col'1
        di as result %-24s "`nm'" _col(30) %14.4f r(sd)
    }
    restore
    matrix drop __boot_acc
    capture matrix drop __brow
end

* ===========================================================================
version 17

mata:
mata set matastrict off

// ---- corrected/uncorrected profile log likelihood at sigma vector s --------
real scalar exactnl_ll(real rowvector s, real matrix M, real scalar N,
                       real rowvector SUM, real scalar jac)
{
    real scalar G, ssr, i, ll
    real colvector m0
    real matrix Mg
    G  = cols(s)
    m0 = M[(2::(G+1)), 1]
    Mg = M[(2::(G+1)), (2::(G+1))]
    ssr = M[1,1] - 2*(s*m0) + (s*Mg*s')
    if (ssr <= 0) {
        return(-1e30)
    }
    ll = -N/2*ln(ssr/N)
    if (jac) {
        for (i=1; i<=G; i++) {
            ll = ll + SUM[i]*ln(1-s[i])
        }
    }
    return(ll)
}

// ---- exhaustive grid for G = 1, 2 (matches the tab4 do-file exactly) -------
real rowvector exactnl_grid(real matrix M, real scalar N, real rowvector SUM,
                            real scalar lo, real scalar hi, real scalar step,
                            real scalar jac, real scalar G)
{
    real scalar ns, i, j, sh, sp, ssr, ll, best
    real scalar a0, aHH, aPP, a0H, a0P, aHP
    real rowvector bs
    ns   = round((hi-lo)/step)
    best = -1e30
    bs   = J(1, G, lo)
    a0   = M[1,1]
    if (G == 1) {
        aHH = M[2,2]
        a0H = M[1,2]
        for (i=0; i<=ns; i++) {
            sh  = lo + i*step
            ssr = a0 - 2*sh*a0H + sh*sh*aHH
            if (ssr > 0) {
                ll = -N/2*ln(ssr/N)
                if (jac) {
                    ll = ll + SUM[1]*ln(1-sh)
                }
                if (ll > best) {
                    best = ll
                    bs   = sh
                }
            }
        }
    }
    else {
        aHH = M[2,2]
        aPP = M[3,3]
        a0H = M[1,2]
        a0P = M[1,3]
        aHP = M[2,3]
        for (i=0; i<=ns; i++) {
            sh = lo + i*step
            for (j=0; j<=ns; j++) {
                sp  = lo + j*step
                ssr = a0 - 2*sh*a0H - 2*sp*a0P + sh*sh*aHH + sp*sp*aPP + 2*sh*sp*aHP
                if (ssr > 0) {
                    ll = -N/2*ln(ssr/N)
                    if (jac) {
                        ll = ll + SUM[1]*ln(1-sh) + SUM[2]*ln(1-sp)
                    }
                    if (ll > best) {
                        best = ll
                        bs   = (sh, sp)
                    }
                }
            }
        }
    }
    return(bs)
}

// ---- coordinate golden-section ascent for G >= 3 ---------------------------
real rowvector exactnl_coord(real matrix M, real scalar N, real rowvector SUM,
                             real scalar lo, real scalar hi, real scalar step,
                             real scalar jac, real scalar G)
{
    real scalar g, it, gr, a, b, c, d, fc, fd, tol, change, oldll, newll
    real rowvector s
    s = J(1, G, 0.2)
    gr = (sqrt(5)-1)/2
    tol = 1e-7
    oldll = exactnl_ll(s, M, N, SUM, jac)
    for (it=1; it<=200; it++) {
        for (g=1; g<=G; g++) {
            a = lo
            b = hi
            c = b - gr*(b-a)
            d = a + gr*(b-a)
            s[g] = c
            fc = exactnl_ll(s, M, N, SUM, jac)
            s[g] = d
            fd = exactnl_ll(s, M, N, SUM, jac)
            while ((b-a) > tol) {
                if (fc > fd) {
                    b = d
                    d = c
                    fd = fc
                    c = b - gr*(b-a)
                    s[g] = c
                    fc = exactnl_ll(s, M, N, SUM, jac)
                }
                else {
                    a = c
                    c = d
                    fc = fd
                    d = a + gr*(b-a)
                    s[g] = d
                    fd = exactnl_ll(s, M, N, SUM, jac)
                }
            }
            s[g] = (a+b)/2
        }
        newll = exactnl_ll(s, M, N, SUM, jac)
        change = newll - oldll
        oldll = newll
        if (change < tol) {
            break
        }
    }
    return(s)
}

// ---- numeric Hessian of the corrected profile LL at s ----------------------
real matrix exactnl_hess(real rowvector s, real matrix M, real scalar N,
                         real rowvector SUM, real scalar lo, real scalar hi)
{
    real scalar G, i, j, h, fpp, fpm, fmp, fmm
    real matrix H
    real rowvector sp
    G = cols(s)
    H = J(G, G, 0)
    for (i=1; i<=G; i++) {
        for (j=i; j<=G; j++) {
            h = 1e-4
            if (1 - s[i] < 5*h) {
                h = (1 - s[i])/5
            }
            if (1 - s[j] < 5*h) {
                h = min((h, (1 - s[j])/5))
            }
            if (s[i] < 5*h) {
                h = min((h, s[i]/5))
            }
            if (s[j] < 5*h) {
                h = min((h, s[j]/5))
            }
            sp = s
            sp[i] = sp[i]+h
            sp[j] = sp[j]+h
            fpp = exactnl_ll(sp, M, N, SUM, 1)
            sp = s
            sp[i] = sp[i]+h
            sp[j] = sp[j]-h
            fpm = exactnl_ll(sp, M, N, SUM, 1)
            sp = s
            sp[i] = sp[i]-h
            sp[j] = sp[j]+h
            fmp = exactnl_ll(sp, M, N, SUM, 1)
            sp = s
            sp[i] = sp[i]-h
            sp[j] = sp[j]-h
            fmm = exactnl_ll(sp, M, N, SUM, 1)
            H[i,j] = (fpp - fpm - fmp + fmm)/(4*h*h)
            H[j,i] = H[i,j]
        }
    }
    return(H)
}

// ---- driver called from the ado --------------------------------------------
void exactnl_fit(string scalar resvars, string scalar touse, string scalar sumname,
                 real scalar N, real scalar lo, real scalar hi, real scalar step,
                 real scalar jac, real scalar G)
{
    real matrix E, M, H, Vs
    real rowvector SUM, sig, se
    real scalar LLbest, LL0
    E   = st_data(., tokens(resvars), touse)
    M   = quadcross(E, E)
    SUM = st_matrix(sumname)
    if (G <= 2) {
        sig = exactnl_grid(M, N, SUM, lo, hi, step, jac, G)
    }
    else {
        sig = exactnl_coord(M, N, SUM, lo, hi, step, jac, G)
    }
    LLbest = exactnl_ll(sig, M, N, SUM, 1)
    LL0    = exactnl_ll(J(1,G,0), M, N, SUM, 1)
    H  = exactnl_hess(sig, M, N, SUM, lo, hi)
    Vs = invsym(-H)
    se = sqrt(diagonal(Vs))'
    st_matrix("__exactnl_sig", sig)
    st_matrix("__exactnl_sigse", se)
    st_numscalar("__exactnl_LLbest", LLbest)
    st_numscalar("__exactnl_LL0", LL0)
}

end
