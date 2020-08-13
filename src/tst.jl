using ForwardDiff
using Plots

#-----------------------------------------------------------------------------
# ---- Cross Sections
#-----------------------------------------------------------------------------

model = KineticMixing(1e2, 1e3, 1.0, 1.0)
cross_section_xx_bb(2 * 6500, model) / (2.56819e-9)


cross_section_xx_dd(4e4, model)
cross_section_xx_ss(4e4, model)
cross_section_xx_bb(4e4, model)
cross_section_xx_ee(4e4, model)
cross_section_xx_μμ(4e4, model)
cross_section_xx_ττ(4e4, model)
cross_section_xx_νν(4e4, model)
cross_section_xx_uu(4e4, model)
cross_section_xx_cc(4e4, model)
cross_section_xx_tt(4e4, model)
cross_section_xx_vv(4e4, model)

@code_typed sigma_xx_tot(4e2, model)

begin
    model = KineticMixing(1e2, 1e3, 1.0, 1e-3)
    cmes = 10 .^ LinRange(log10(2.001model.mx), log10(100model.mx), 1000)
    cs = [cross_section_xx_tot(cme, model) for cme in cmes]

    plot(cmes, cs, xaxis = :log, yaxis = :log)
end

#-----------------------------------------------------------------------------
# ---- Thermal Cross Sections ------------------------------------------------
#-----------------------------------------------------------------------------

begin
    model = KineticMixing(0.9e3, 1e3, 1.0, 1e-3)
    xs = 10 .^ LinRange(-1.0, 3.0, 150)
    tcs = [thermal_cross_section_xx_tot(x, model) for x in xs]

    plot(xs, tcs, xaxis = :log, yaxis = :log)
end

#-----------------------------------------------------------------------------
# ---- other ------------------------------------------------
#-----------------------------------------------------------------------------

begin
    Ts = 10 .^ LinRange(-5.0, 5.0, 200)
    hs = sm_heff.(Ts)
    gs = sm_geff.(Ts)
    gstars = sm_sqrt_gstar.(Ts)

    plot(Ts, hs, label = "heff")
    plot!(Ts, gs, label = "geff")
    plot!(Ts, gstars, label = "√g*")
    xaxis!(:log)
end


#-----------------------------------------------------------------------------
# ---- other ------------------------------------------------
#-----------------------------------------------------------------------------

using ODEInterfaceDiffEq

@code_warntype sm_geff(1.0)

xspan = (0.1, 1000.0)
logxspan = (log(xspan[1]), log(xspan[2]))
u0 = [weq(model.mx / xspan[1], model.mx, 2.0, 1)]
f = ODEFunction(boltzmann!, jac = boltzmann_jac!)

prob = ODEProblem(boltzmann!, u0, logxspan, model)

sol_radau = solve(prob, alg = radau5(), reltol = 1e-9, abstol = 1e-9)
sol_rodas = solve(prob, alg = rodas(), reltol = 1e-9, abstol = 1e-9)
sol_seulex = solve(prob, alg = seulex(), reltol = 1e-9, abstol = 1e-9)
sol_ros23 = solve(prob, alg = Rosenbrock23(), reltol = 1e-9, abstol = 1e-9)


xspan = (BigFloat(0.1), BigFloat(1000.0))
logxspan = (log(xspan[1]), log(xspan[2]))
u0 = [BigFloat(weq(model.mx / Float64(xspan[1]), model.mx, 2.0, 1))]
f = ODEFunction(boltzmann!, jac = boltzmann_jac!)
sol_radauIIA = solve(prob, alg = RadauIIA5(), reltol = BigFloat(1e-20), abstol = 1e-20)


plot(sol_radau)
plot!(sol_rodas)
plot!(sol_seulex)
plot!(sol_ros23)
plot!(sol_radauIIA)

sol_radauIIA.u[end]
