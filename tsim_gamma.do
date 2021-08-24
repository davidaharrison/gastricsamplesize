
capture program drop tsim_gamma
program tsim_gamma, rclass
	args a b d ninf reps ss
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
		*generate ventilation days, using gamma distribution
		gen days=rgamma(`a', `b')
		replace days=30 if days>30
		replace days=31 if death==1
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