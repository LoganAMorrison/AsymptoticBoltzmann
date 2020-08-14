using ForwardDiff
using Plots

#-----------------------------------------------------------------------------
# ---- Cross Sections
#-----------------------------------------------------------------------------

model = KineticMixing(1e2, 1e3, 1.0, 1e-3)

sigma_xx_dd(4e4, model)
sigma_xx_ss(4e4, model)
sigma_xx_bb(4e4, model)
sigma_xx_ee(4e4, model)
sigma_xx_μμ(4e4, model)
sigma_xx_ττ(4e4, model)
sigma_xx_νν(4e4, model)
sigma_xx_uu(4e4, model)
sigma_xx_cc(4e4, model)
sigma_xx_tt(4e4, model)
sigma_xx_vv(4e4, model)

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

MV = 1e3
EPS = 1e-3
GVXX = 1.0
rs = LinRange(0.1, 0.6, 50)
mxs = MV .* rs

rds = Array{Tuple{Float64, Float64, Float64}}(undef, length(rs))



for i in 1:length(mxs)
    mx = mxs[i]
    model = KineticMixing(mx, MV, GVXX, EPS)
    α = (MV - 2*mx) / mx
    try
    rds[i] = (
        compute_relic_density_mpu(model),
        compute_relic_density_mpu(model, α),
        compute_relic_density_radau(model)
    )
catch
    rds[i] = (
        NaN,NaN,NaN
    )
end
end



#-----------------------------------------------------------------------------
# ---- other ------------------------------------------------
#-----------------------------------------------------------------------------

begin
    model = KineticMixing(3e2, 1e3, 1.0, 1e-3)

    xs = LinRange(0.1, 50.0, 500)
    fs = [compute_f(x, model) for x in xs]
    qs = [compute_q(x, model) for x in xs]
    ps = [compute_p(x, model) for x in xs]

    plot(xs, qs, label = "λQ(x)")
    plot!(xs, abs.(ps), label = "P(x)")
    plot!(xs, fs * 8e-9, label = "λf(x)")
    xlabel!("x")
    yaxis!(:log)
end
