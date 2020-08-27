include("../src/AsymptoticBoltzmann.jl")


function mx_scan(mxs, mv, gvxx, eps)
    rds = Array{Float64,2}(undef, length(mxs), 3)

    for i = 1:length(mxs)
        mx = mxs[i]
        model = AsymptoticBoltzmann.KineticMixing(mx, mv, gvxx, eps)

        try
            rds[i, 1] = AsymptoticBoltzmann.relic_density(
                model,
                AsymptoticBoltzmann.ODE(reltol = 1e-12),
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
            α = max((mv - 2mx) / mx, 0.0)
            rds[i, 3] = AsymptoticBoltzmann.relic_density(
                model,
                AsymptoticBoltzmann.MPU(α=α),
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
    rs = LinRange(0.1, 0.6, 150)
    mxs = rs .* MV

    mx_scan(mxs, MV, GVXX, EPS)
end

using Plots


plot(LinRange(0.1, 0.6, 150),
    abs.(rds[:,1] - rds[:,2]) ./ rds[:,1],
    yaxis=:log
)
plot!(LinRange(0.1, 0.6, 150),
    abs.(rds[:,1] - rds[:,3]) ./ rds[:,1],
    yaxis=:log
)
