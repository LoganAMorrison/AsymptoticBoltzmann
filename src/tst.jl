using ForwardDiff
using Plots

#-----------------------------------------------------------------------------
# ---- Cross Sections
#-----------------------------------------------------------------------------

model = KineticMixing(1e2, 1e3, 1.0, 1e-3)

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

@code_warntype dm_annihilation_cross_section(4e2, model)

model.width_v

begin
    model = KineticMixing(350.0, 1e3, 1.0, 1e-2)
    cmes = 10 .^ LinRange(log10(2.001model.mx), log10(4000), 1000)
    cs = [dm_annihilation_cross_section(cme, model) for cme in cmes]

    plot(cmes, cs, xaxis = :log, yaxis = :log)
end

#-----------------------------------------------------------------------------
# ---- Thermal Cross Sections ------------------------------------------------
#-----------------------------------------------------------------------------

begin
    model = KineticMixing(9e2, 1e3, 1.0, 1e-3)
    xs = 10 .^ LinRange(-1.0, 3.0, 500)
    tcs = [dm_thermal_cross_section(x, model) for x in xs]

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

model = KineticMixing(1e2, 1e3, 1.0, 1e-3)

relic_density(model, MPU())
relic_density(model, ODE(reltol = 1e-12))


#-----------------------------------------------------------------------------
# ---- other ------------------------------------------------
#-----------------------------------------------------------------------------

begin
    model = KineticMixing(3e2, 1e3, 1.0, 1e-3)

    xs = LinRange(1e-3, 50.0, 1000)
    fs = [compute_f(x, model) for x in xs]
    qs = [compute_q(x, model) for x in xs]
    ps = [compute_p(x, model) for x in xs]

    plot(xs, qs, label = "λQ(x)")
    plot!(xs, abs.(ps), label = "|P(x)|")
    plot!(xs, fs * 8e-9, label = "λf(x)")
    xlabel!("x")
    yaxis!(:log)
end
