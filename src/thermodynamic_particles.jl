function neq(T::Real, mass::Real, g::Real, spin2::Int)
    x = mass / T
    eta = spin2 % 2 == 0 ? 1 : -1
    g * x^2 * T^3 / (2pi^2) * besselk(2, x)
end

function yeq(T::Real, mass::Real, g::Real, spin2::Int)
    x = mass / T
    eta = spin2 % 2 == 0 ? 1 : -1

    45g * x^2 / (4π^4 * sm_heff(T)) * besselk(2, x)
end

function weq(T::Real, mass::Real, g::Real, spin2::Int)
    x = mass / T
    eta = spin2 % 2 == 0 ? 1 : -1

    log(45g * x^2 / (4π^4 * sm_heff(T)) * besselk(2, x))
end

function weq(T::Float64, mass::Real, g::Real, spin2::Int)
    x = mass / T
    eta = spin2 % 2 == 0 ? 1 : -1

    -x + log(45g * x^2 / (4π^4 * sm_heff(T)) * besselkx(2, x))
end
