
capture program drop roc_lognormal
program roc_lognormal, rclass
	args mC sC dC mI sI dI
preserve
quietly {					
    clear
		
	set obs 31
	gen x=_n
	*probability of x in control arm
	gen double pr_x_C=`dC' if x==31
	replace pr_x_C=(1-`dC')*(normal((ln(x)-`mC')/`sC')) if x==1
	replace pr_x_C=(1-`dC')*(normal((ln(x)-`mC')/`sC')-normal((ln(x-1)-`mC')/`sC')) if x>1 & x<30
	replace pr_x_C=(1-`dC')*(1-normal((ln(29)-`mC')/`sC')) if x==30

/*	*cumulative probability of x in control arm
	gen double cpr_x_C=1 if x==31
	replace cpr_x_C=1-`dC' if x==30
	replace cpr_x_C=(1-`dC')*(normal((ln(x)-`mC')/`sC')) if x==1
	replace cpr_x_C=(1-`dC')*(normal((ln(x)-`mC')/`sC')) if x>1 & x<30

	*probability of x in intervention arm
	gen double pr_x_I=`dI' if x==31
	replace pr_x_I=(1-`dI')*(normal((ln(x)-`mI')/`sI')) if x==1
	replace pr_x_I=(1-`dI')*(normal((ln(x)-`mI')/`sI')-normal((ln(x-1)-`mI')/`sI')) if x>1 & x<30
	replace pr_x_I=(1-`dI')*(1-normal((ln(29)-`mI')/`sI')) if x==30
*/
	*cumulative probability of x in intervention arm
	gen double cpr_x_I=1 if x==31
	replace cpr_x_I=1-`dI' if x==30
	replace cpr_x_I=(1-`dI')*(normal((ln(x)-`mI')/`sI')) if x==1
	replace cpr_x_I=(1-`dI')*(normal((ln(x)-`mI')/`sI')) if x>1 & x<30
	
	gen double p=pr_x_C*cpr_x_I/2 if x==1
	replace p=pr_x_C*(cpr_x_I + cpr_x_I [_n-1])/2 if x>1
	total p
	local auc=el(r(table),1,1)

}
di `auc'
return scalar auc=`auc'
restore
end