function compute_f(x::Real, model)
    mx = dark_matter_mass(model)
    T = mx / x
    sv = thermal_cross_section(x, model)
    sqrt(π / 45) * mx * M_PLANK / x^2 * sm_sqrt_gstar(T) * sv
end

function compute_q(x::Real, model)
    mx = dark_matter_mass(model)
    T = mx / x
    sv = thermal_cross_section(x, model)
    ye = yeq(T, model)
    sqrt(π / 45) * mx * M_PLANK / x^2 * sm_sqrt_gstar(T) * sv * ye
end

function compute_dq(x::Real, model)
    first(ForwardDiff.partials(compute_f(ForwardDiff.Dual{Float64}(x, 1.0), model)))
end

function compute_p(x::Real, model)
    mx = dark_matter_mass(model)
    T = mx / x
    sv = thermal_cross_section(x, model)
    f(xx) = (
        sqrt(π / 45) * mx * M_PLANK / xx^2 *
        sm_sqrt_gstar(mx / xx) *
        thermal_cross_section(xx, model)
    )

    val = f(ForwardDiff.Dual{ForwardDiff.Dual{Float64}}(
        ForwardDiff.Dual{Float64}(x, 1.0),
        ForwardDiff.Dual{Float64}(1.0, 0.0),
    ))

    f = ForwardDiff.value(ForwardDiff.value(val))
    df = first(ForwardDiff.partials(ForwardDiff.value(val)))
    d2f = first(ForwardDiff.partials(first(ForwardDiff.partials(val))))

    3 * (df / f)^2 / 4 - (d2f / f) / 2
end

function compute_xf(model)
    f(xf) = compute_q(xf, model) - 0.41309880210998464
    fzero(f, 0.1, 100.0)
end

function compute_δx(α, model)
    α == zero(typeof(α)) && return one(typeof(α))
    (log(gamma((2 + 3α) / (2 + 2α)) / gamma((2 + α) / (2 + 2α))) / α +
    (1.5772156649015329 + log(1 + α)) / (1 + α))
end

function compute_relic_density_mpu(model, α=0.0)
    mx = dark_matter_mass(model)
    xf = compute_xf(model)
    δx = compute_δx(α, model)
    function f(x)
        try
            return compute_f(x, model)
        catch
            return 0.0
        end
    end
    mx * S_TODAY / RHO_CRIT / quadgk(f, xf - δx, Inf)[1]
end
