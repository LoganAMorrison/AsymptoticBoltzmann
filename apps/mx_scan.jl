include("../src/AsymptoticBoltzmann.jl")
using Plots

function mx_scan(mxs, mv, gvxx, eps)
    rds = Array{Float64,2}(undef, length(mxs), 3)

    Threads.@threads for i = 1:length(mxs)
        mx = mxs[i]
        model = AsymptoticBoltzmann.KineticMixing(mx, mv, gvxx, eps)

        try
            rds[i, 1] = AsymptoticBoltzmann.relic_density(
                model,
                AsymptoticBoltzmann.ODE(
                    xstart = 1,
                    xend = 5e4,
                    reltol = 1e-14,
                    abstol = 2e-19,
                ),
            )
        catch
            rds[i, 1] = NaN
        end

        try
            rds[i, 2] = AsymptoticBoltzmann.relic_density(
                model,
                AsymptoticBoltzmann.MPU(),
            )
        catch
            rds[i, 2] = NaN
        end

        try
            if 2mx < mv
                α = max((mv - 2mx) / mx, 0.0)
            else
                α = 2 * (mv - mx) / mx
            end
            rds[i, 3] = AsymptoticBoltzmann.relic_density(
                model,
                AsymptoticBoltzmann.MPU(α = max(α, 0.0)),
            )
        catch
            rds[i, 3] = NaN
        end
    end
    rds
end

rds = begin
    MV = 1e3
    EPS = 1e-3
    GVXX = 1.0
    rs = 10 .^ LinRange(log10(0.1), log10(1.2), 500)
    mxs = rs .* MV

    mx_scan(mxs, MV, GVXX, EPS)
end

xfs = begin
    MV = 1e3
    EPS = 1e-3
    GVXX = 1.0
    rs = 10 .^ LinRange(log10(0.1), log10(1.2), 500)
    mxs = rs .* MV

    _xfs = zeros(eltype(mxs), length(mxs))

    for (i, mx) in enumerate(mxs)
        model = AsymptoticBoltzmann.KineticMixing(mx, MV, GVXX, EPS)
        _xfs[i] = AsymptoticBoltzmann.compute_xf(model)
    end
    _xfs
end

tcs = begin

    MV = 1e3
    EPS = 1e-3
    GVXX = 1.0
    MX = MV / 2 * 1.1
    model = AsymptoticBoltzmann.KineticMixing(MX, MV, GVXX, EPS)

    xs = 10 .^ LinRange(-1.0, 2.0, 500)
    _tcs = zeros(eltype(mxs), length(mxs))
    for (i,x) in enumerate(xs)
        _tcs[i] = AsymptoticBoltzmann.dm_thermal_cross_section(x,model)
    end
    _tcs
end

cs = begin

    MV = 1e3
    EPS = 1e-3
    GVXX = 1.0
    MX = 100.0
    model = AsymptoticBoltzmann.KineticMixing(MX, MV, GVXX, EPS)

    Qs = 10 .^ LinRange(log10(2.001*MX), log10(50*MX), 500)
    _cs = zeros(eltype(mxs), length(mxs))
    for (i,Q) in enumerate(Qs)
        _cs[i] = AsymptoticBoltzmann.dm_annihilation_cross_section(Q,model)
    end
    _cs
end

log.(rds[:, 1])


plot(rs, abs.(rds[:, 1] - rds[:, 2]) ./ rds[:, 1], yaxis = :log)
plot!(rs, abs.(rds[:, 1] - rds[:, 3]) ./ rds[:, 1], yaxis = :log)

plot(rs, log.(rds[:, 1]))

@show log.(rds[:, 1])


@show cs
