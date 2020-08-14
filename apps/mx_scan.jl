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
            rds[i, 3] = AsymptoticBoltzmann.relic_density(
                model,
                AsymptoticBoltzmann.MPU(α = (mv - 2mx) / mx),
            )
        catch
            rds[i, 3] = NaN
        end
    end
    rds
end

begin
    MV = 1e3
    EPS = 1e-3
    GVXX = 1.0
    rs = LinRange(0.1, 0.6, 20)
    mxs = rs .* MV

    mx_scan(mxs, MV, GVXX, EPS)
end
