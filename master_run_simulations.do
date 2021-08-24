/* Project: GASTRIC*/ 
/* Purpose of do-file: create program to simulate trial results, using log-normal dist */ 
/* Author: Karen Thomas (copied from DH's simulations for Big Oxy)*/ 
/* Date: march 2019*/ 
/* Required input: */
/* Output: log file*/

version 16.1
set more off
cap log close
set seed 79512

clear

cd "S:\CTU\14. Study proposals\GASTRIC\Sample size\Jan2021"

do "roc_lognormal.do"
do "tsim_lognormal.do"
do "roc_gamma.do"
do "tsim_gamma.do"

*****explore log-normal - median 2.9, lower quartile 1.0*****
*************************************************************
clear
gen str50 control_dist_daysv=""
gen double control_mort=.
gen str50 int_dist_daysv=""
gen double int_mort=.
gen double rho=.
gen int ss90=.


local medC=2.9
local lqC=1.0

local mC=ln(`medC')
local sC=(ln(`lqC')-ln(`medC'))/(invnormal(0.25))
local dC=0.042

forvalues mh=9(1)10 {
	local mdiff=`mh'/24
	local qdiff=(`lqC'*(`medC'+`mdiff')/`medC')-`lqC'

	local mI=ln(`medC'+`mdiff')
	local sI=(ln(`lqC'+`qdiff')-ln(`medC'+`mdiff'))/(invnormal(0.25))
	local new=_N+1	
	set obs `new'
	local dI=0.042
	replace control_dist_daysv="Log-normal, median `medC' lower quartile `lqC'" if _n==_N
	replace control_mort=4.2 if _n==_N
	replace int_dist_daysv="Log-normal, median `=round(`medC'+`mdiff',0.01)' lower quartile `=round(`lqC'+`qdiff',0.01)'" if _n==_N
	replace int_mort=4.2 if _n==_N
	roc_lognormal `mC' `sC' `dC' `mI' `sI' `dI'
	replace rho=r(auc) if _n==_N
	local r=round(r(auc),0.001)
	local p=0
	local s=0
	forvalues ss=4000(100)5800 {
		if `p'<90 {
				tsim_lognormal `mC' `sC' `dC' `r' 1000 `ss'
				local p=r(pwr)
				if `p'>=90 {
					local s=`ss'
				}
			}
		}
	replace ss90=`s' if _n==_N
	append using simresults
	save simresults, replace
	}

*****explore log-normal - median 3.0, lower quartile 2.0*****
*************************************************************
local medC=3.0
local lqC=2.0

local mC=ln(`medC')
local sC=(ln(`lqC')-ln(`medC'))/(invnormal(0.25))
local dC=0.042

forvalues mh=4(0.5)6 {
	local mdiff=`mh'/24
	local qdiff=(`lqC'*(`medC'+`mdiff')/`medC')-`lqC'

	local mI=ln(`medC'+`mdiff')
	local sI=(ln(`lqC'+`qdiff')-ln(`medC'+`mdiff'))/(invnormal(0.25))
	local new=_N+1	
	set obs `new'
	local dI=0.042
	replace control_dist_daysv="Log-normal, median `medC' lower quartile `lqC'" if _n==_N
	replace control_mort=4.2 if _n==_N
	replace int_dist_daysv="Log-normal, median `=round(`medC'+`mdiff',0.01)' lower quartile `=round(`lqC'+`qdiff',0.01)'" if _n==_N
	replace int_mort=4.2 if _n==_N
	roc_lognormal `mC' `sC' `dC' `mI' `sI' `dI'
	replace rho=r(auc) if _n==_N
	local r=round(r(auc),0.001)
	local p=0
	local s=0
	forvalues ss=2600(100)4800 {
		if `p'<90 {
				tsim_lognormal `mC' `sC' `dC' `r' 10000 `ss'
				local p=r(pwr)
				if `p'>=90 {
					local s=`ss'
				}
			}
		}
	replace ss90=`s' if _n==_N
	}

		
format rho %9.3f
log using "simtable.txt", text replace
list, noobs sep(0)
log close

*****explore log-normal - median 2.9, lower quartile 1.0, with increasing mortality*****
*************************************************************
local medC=2.9
local lqC=1.0

local mC=ln(`medC')
local sC=(ln(`lqC')-ln(`medC'))/(invnormal(0.25))

local mdiff=0.5
local qdiff=(`lqC'*(`medC'+`mdiff')/`medC')-`lqC'

local mI=ln(`medC'+`mdiff')
local sI=(ln(`lqC'+`qdiff')-ln(`medC'+`mdiff'))/(invnormal(0.25))

forvalues d=0.042(0.01)0.092 {
	local dC=`d'
	local dI=`d'
	di "Mortality `=100*`d''"
	roc_lognormal `mC' `sC' `dC' `mI' `sI' `dI'
	local r=round(r(auc),0.001)
	*tsim_lognormal `mC' `sC' `dC' `r' 1000 4000
}


	

	
*look at differences in mortality

local mC=ln(2.9)
local sC=(ln(1.0 )-ln(2.9))/(invnormal(0.25))
local dC=0.042
log using "difference in mortality.txt", text replace
forvalues diff=0.005(0.005)0.05 {
	di "Control arm `dC' Intervention arm + `=`diff'''"
	local mI=ln(2.9)
	local sI=(ln(1.0)-ln(2.9))/(invnormal(0.25))
	local dI=`dC' + `diff'

	roc_lognormal `mC' `sC' `dC' `mI' `sI' `dI'

	}
log close

*differences in energy targets
power twomeans 0.7 0.75, sd1(0.2 0.25 0.3 0.35 0.4) sd2(0.2 0.25 0.3 0.35 0.4) power(0.9)


*now gamma (a is shape, b is scale - the mean is a/b so if we keep b fixed at 1, increase a in increments)
local aC=2.9
local bC=1.0
local dC=0.042



forvalues diff=0.5(0.1)1.5 {
	local new=_N+1	
	set obs `new'
	local aI=2.9+`diff'
	local bI=1.0
	local dI=0.042

	replace control_dist_daysv="Gamma(2.9, 1.0)" if _n==_N
	replace control_mort=4.2 if _n==_N
	replace int_dist_daysv="Gamma(`=2.9+`diff'', 1.0)" if _n==_N
	replace int_mort=4.2 if _n==_N
	roc_gamma `aC' `bC' `dC' `aI' `bI' `dI'
	replace rho=r(auc) if _n==_N

	local r=round(r(auc),0.01)
	local p=0
	local s=0
	forvalues ss=1000(200)3000 {
		if `p'<90 {
				tsim_gamma `aC' `bC' `dC' `r' 1000 `ss'
				local p=r(pwr)
				if `p'>=90 {
					local s=`ss'
				}
			}
		}
	replace ss90=`s' if _n==_N
	}
save as results, replace

exit




