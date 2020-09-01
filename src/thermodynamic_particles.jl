function neq(T::Real, mass::Real, g::Real, spin2::Int)
    x = mass / T
    η = spin2 % 2 == 0 ? 1 : -1
    # g * x^2 * T^3 / (2pi^2) * besselk(2, x)

    f(y) = y * sqrt(y^2 - 1) / (exp(y * x) - η)
    g * mass^3 / (2π^2) * quadgk(f, 1, Inf)[1]
end

function yeq(T::Real, mass::Real, g::Real, spin2::Int)
    x = mass / T
    η = spin2 % 2 == 0 ? 1 : -1
    # 45g * x^2 / (4π^4 * sm_heff(T)) * besselk(2, x)

    f(y) = y * sqrt(y^2 - 1) / (exp(y * x) - η)
    45g * x^3 / (4π^4 * sm_heff(T)) * quadgk(f, 1, Inf)[1]
end

function weq(T::Real, mass::Real, g::Real, spin2::Int)
    x = mass / T
    η = spin2 % 2 == 0 ? 1 : -1

    #log(45g * x^2 / (4π^4 * sm_heff(T)) * besselk(2, x))

    f(y) = y * sqrt(y^2 - 1) * exp(-x * (y - 1)) / (1 - η * exp(-y * x))
    -x + log(45g * x^3 / (4π^4 * sm_heff(T)) * quadgk(f, 1, Inf)[1])
end

using QuadGK
