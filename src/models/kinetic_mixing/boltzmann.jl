function boltzmann!(dw, w, model::KineticMixing, logx)
    x = exp(logx)
    T = model.mx / x
    s = sm_entropy_density(T)
    we = weq(T, model.mx, 2.0, 1)
    pf = -sqrt(π / 45) * M_PLANK * sm_sqrt_gstar(T) * T
    sv = thermal_cross_section_xx_tot(x, model)

    dw[1] = pf * sv * (exp(w[1]) - exp(2we - w[1]))
end

function boltzmann_jac!(J, w, model::KineticMixing, logx)
    x = exp(logx)
    T = model.mx / x
    s = sm_entropy_density(T)
    we = weq(T, model.mx, 2.0, 1)
    pf = -sqrt(π / 45) * M_PLANK * sm_sqrt_gstar(T) * T
    sv = thermal_cross_section_xx_tot(x, model)

    J[1, 1] = pf * sv * (exp(w[1]) + exp(2we - w[1]))
end
