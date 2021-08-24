
capture program drop roc_gamma
program roc_gamma, rclass
	args aC bC dC aI bI dI
preserve
quietly {					
    clear
	
	set obs 31
	gen x=_n
	*probability of x in control arm
	gen double pr_x_C=`dC' if x==31
	replace pr_x_C=(1-`dC')*(gammap(`aC', x/`bC')) if x==1
	replace pr_x_C=(1-`dC')*(gammap(`aC', x/`bC')-gammap(`aC', (x-1)/`bC')) if x>1 & x<30
	replace pr_x_C=(1-`dC')*(1-gammap(`aC', 29/`bC')) if x==30

	*cumulative probability of x in intervention arm
	gen double cpr_x_I=1 if x==31
	replace cpr_x_I=1-`dI' if x==30
	replace cpr_x_I=(1-`dI')*(gammap(`aI', x/`bI')) if x<30
	
	gen double p=pr_x_C*cpr_x_I/2 if x==1
	replace p=pr_x_C*(cpr_x_I + cpr_x_I [_n-1])/2 if x>1
	total p
	local auc=el(r(table),1,1)

}
di `auc'
return scalar auc=`auc'
restore
end