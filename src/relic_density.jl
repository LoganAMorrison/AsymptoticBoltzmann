abstract type AbstractRelicDensitySolver end

# ==========================================================================
# ==== ODE Solving with Radau ==============================================
# ==========================================================================

struct ODE
    alg::OrdinaryDiffEqAlgorithm
    xstart::Float64
    xend::Float64
    reltol::Float64
    abstol::Float64
end

function ODE(;
    alg = RadauIIA5(),
    xstart = 1,
    xend = 1e4,
    reltol = 1e-20,
    abstol = 1e-20,
)
    ODE(alg, xstart, xend, reltol, abstol)
end

function boltzmann!(dw, w, model::AbstractDMModel, logx)
    mx = dm_mass(model)
    x = exp(logx)
    T = mx / x
    s = sm_entropy_density(T)
    we = weq(T, mx, dm_dof(model), dm_spin2(model))
    pf = -sqrt(π / 45) * M_PLANK * sm_sqrt_gstar(T) * T
    sv = dm_thermal_cross_section(x, model)

    dw[1] = pf * sv * (exp(w[1]) - exp(2we - w[1]))
end

function boltzmann_jac!(J, w, model::AbstractDMModel, logx)
    mx = dm_mass(model)
    x = exp(logx)
    T = mx / x
    s = sm_entropy_density(T)
    we = weq(T, mx, dm_dof(model), dm_spin2(model))
    pf = -sqrt(π / 45) * M_PLANK * sm_sqrt_gstar(T) * T
    sv = dm_thermal_cross_section(x, model)

    J[1, 1] = pf * sv * (exp(w[1]) + exp(2we - w[1]))
end

function relic_density(model::AbstractDMModel, ode::ODE)
    mx = dm_mass(model)
    logxspan = (log(ode.xstart), log(ode.xend))
    u0 = [weq(mx / ode.xstart, mx, 2.0, 1)]
    f = ODEFunction(boltzmann!, jac = boltzmann_jac!)
    prob = ODEProblem(f, u0, logxspan, model)
    sol = solve(prob, alg = ode.alg, reltol = ode.reltol, abstol = ode.abstol)
    exp(sol.u[end][1]) * mx * S_TODAY / RHO_CRIT
end

# ==========================================================================
# ==== Morrison-Patel-Ulbricht =============================================
# ==========================================================================

struct MPU <: AbstractRelicDensitySolver
    α::Float64
end

function MPU(; α = 0.0)
    MPU(α)
end

function compute_f(x::Real, model::AbstractDMModel)
    mx = dm_mass(model)
    T = mx / x

    sqrt(π / 45) * mx * M_PLANK / x^2 *
    sm_sqrt_gstar(T) *
    dm_thermal_cross_section(x, model)
end

function compute_q(x::Real, model::AbstractDMModel)
    mx = dm_mass(model)
    T = mx / x
    sv = dm_thermal_cross_section(x, model)
    ye = yeq(T, mx, dm_dof(model), dm_spin2(model))
    sqrt(π / 45) * mx * M_PLANK / x^2 * sm_sqrt_gstar(T) * sv * ye
end

function compute_dq(x::Real, model::AbstractDMModel)
    first(ForwardDiff.partials(compute_f(
        ForwardDiff.Dual{Float64}(x, 1.0),
        model,
    )))
end

function compute_p(x::Real, model::AbstractDMModel)
    mx = dm_mass(model)
    T = mx / x
    sv = dm_thermal_cross_section(x, model)
    f(xx) = compute_f(xx, model)

    val = f(ForwardDiff.Dual{ForwardDiff.Dual{Float64}}(
        ForwardDiff.Dual{Float64}(x, 1.0),
        ForwardDiff.Dual{Float64}(1.0, 0.0),
    ))

    f = ForwardDiff.value(ForwardDiff.value(val))
    df = first(ForwardDiff.partials(ForwardDiff.value(val)))
    d2f = first(ForwardDiff.partials(first(ForwardDiff.partials(val))))

    3 * (df / f)^2 / 4 - (d2f / f) / 2
end

function compute_xf(model::AbstractDMModel)
    f(xf) = compute_q(xf, model) - 0.41309880210998464
    fzero(f, 0.1, 100.0)
end

function compute_δx(α, model::AbstractDMModel)
    α == zero(typeof(α)) && return one(typeof(α))
    (
        log(gamma((2 + 3α) / (2 + 2α)) / gamma((2 + α) / (2 + 2α))) / α +
        (1.5772156649015329 + log(1 + α)) / (1 + α)
    )
end

function relic_density(model::AbstractDMModel, mpu::MPU)
    mx = dm_mass(model)
    xf = compute_xf(model)
    δx = compute_δx(mpu.α, model)
    function f(x)
        try
            return compute_f(x, model)
        catch
            return 0.0
        end
    end
    mx * S_TODAY / RHO_CRIT / quadgk(f, xf - δx, Inf)[1]
end
