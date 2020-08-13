_sm_data = CSV.File(joinpath(@__DIR__, "data.csv"))

_sm_heff = LinearInterpolation(
    log10.(_sm_data.T),
    _sm_data.heff,
    extrapolation_bc = Flat(),
);

_sm_geff = LinearInterpolation(
    log10.(_sm_data.T),
    _sm_data.geff,
    extrapolation_bc = Flat(),
);

_sm_sqrt_gstar = LinearInterpolation(
    log10.(_sm_data.T),
    _sm_data.gstar,
    extrapolation_bc = Flat(),
);

sm_heff(T::Real) = _sm_heff(log10(T))
sm_geff(T::Real) = _sm_geff(log10(T))
sm_sqrt_gstar(T::Real) = _sm_sqrt_gstar(log10(T))

sm_entropy_density(T::Real) = 2π^2 / 45 * sm_heff(T) * T^3
sm_energy_density(T::Real) = π^2 / 30 * sm_geff(T) * T^4
