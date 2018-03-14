program define ch_weak, rclass
	syntax [if] [in], p(numlist) beta_range(numlist)   y(name) x(name) z(varlist) [controls(varlist)] [absorb(varname)] [WEIGHT_var(string)] [cluster(varlist)]
	local x `x'
	local y `y'
    if "`weight_var'" != "" {
        local weight "[aw=`weight_var']"
    }
    if "`cluster'" != "" {
        local se "cluster(`cluster')"
    }
    else {
        local se "robust"
    }

    qui areg `y' `controls' `weight', absorb(`absorb')
    tempvar y_res
    predict `y_res', res
    qui areg `x' `controls' `weight', absorb(`absorb')
    tempvar x_res
    predict `x_res', res
    
    local beta_min = .
    local beta_max = .
    foreach beta in `beta_range' {
        tempvar y_tilde
        qui gen `y_tilde' = `y_res' - `x_res'*`beta'
        qui reg `y_tilde' `z' `weight', `se'
        qui test `z' 
        if `beta_min' == . & r(p) > `p' {
            local beta_min = `beta'
        }
        if (r(p) > `p'  & (`beta_max' < `beta' | `beta_max' == .)) {
            local beta_max = `beta'            
        }
        mat val = (`beta' \ r(p))
        mat p = (nullmat(p), val)
    }
    return matrix p = p
    return scalar beta_max = `beta_max'
    return scalar beta_min = `beta_min'
end
