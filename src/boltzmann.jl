function boltzmann!(dw, w, model, logx)
    x = exp(logx)
    T = model.mx / x
    s = sm_entropy_density(T)
    we = weq(T, model.mx, 2.0, 1)
    pf = -sqrt(π / 45) * M_PLANK * sm_sqrt_gstar(T) * T
    sv = thermal_cross_section(x, model)

    dw[1] = pf * sv * (exp(w[1]) - exp(2we - w[1]))
end

function boltzmann_jac!(J, w, model, logx)
    x = exp(logx)
    T = model.mx / x
    s = sm_entropy_density(T)
    we = weq(T, model.mx, 2.0, 1)
    pf = -sqrt(π / 45) * M_PLANK * sm_sqrt_gstar(T) * T
    sv = thermal_cross_section(x, model)

    J[1, 1] = pf * sv * (exp(w[1]) + exp(2we - w[1]))
end

function compute_relic_density_radau(model; xspan=(1,1e4), reltol=1e-20, abstol=1e-20)
    mx = dark_matter_mass(model)
    logxspan = (log(xspan[1]), log(xspan[2]))
    u0 = [weq(mx / xspan[1], mx, 2.0, 1)]
    f = ODEFunction(boltzmann!, jac = boltzmann_jac!)
    prob = ODEProblem(boltzmann!, u0, logxspan, model)
    sol_radauIIA = solve(prob, alg = RadauIIA5(), reltol = reltol, abstol = abstol)
    exp(sol_radauIIA.u[end][1]) * mx * S_TODAY / RHO_CRIT
end
