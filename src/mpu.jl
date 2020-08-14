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

end

function compute_yinf(model)

end
