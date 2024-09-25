function compute_time_step!(p::HeatProblem, Δt) 
    (; base, K, M, uⁿ, uⁿ⁺¹) = p 
    (; ch, g, J) = base
    # .nzval assures that structural zeros are NOT dropped (-> needed to apply constraints)

    g .= (M ./ Δt .- K ./ 2) * uⁿ
    J.nzval .= (M.nzval ./ Δt + K.nzval ./ 2)
    apply!(J, g, ch) 
    uⁿ⁺¹ .= J \ g 
    return p
end