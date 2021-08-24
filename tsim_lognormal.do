
capture program drop tsim_lognormal
program tsim_lognormal, rclass
	args m s d ninf reps ss
quietly {					
preserve
local nsuccess=0
forvalues r=1/`reps' {
		clear
		set obs `ss'
		*randomise patients 1:1
		gen byte tx=(_n<(`ss'/2))
		*generate deaths
		gen byte death=(uniform()<`d')
		*generate ventilation days, using log-normal distribution
		gen days=exp(`m'+`s'*rnormal(0,1))
		*round to nearest hour
		replace days=round(24*days,1)
		replace days=30*24 if days>30*24
		replace days=31*24 if death==1
		*test for non-inferiority
		roctab tx days, binomial level(90)
		local tsuccess=(r(lb)>=`ninf')
		local nsuccess=`nsuccess' + `tsuccess'
	}
restore
}
local pwr=100*`nsuccess'/`reps'
di string(`pwr' , "%9.0f")
return scalar pwr=`pwr'
end