_sm_data = CSV.File(joinpath(@__DIR__, "smdof.csv"))

_sm_heff = CubicSplineInterpolation(
    LinRange(-4.7, 4.1, 276),
    _sm_data.Heff;
    extrapolation_bc = Flat(),
);

_sm_geff = CubicSplineInterpolation(
    LinRange(-4.7, 4.1, 276),
    _sm_data.Geff,
    extrapolation_bc = Flat(),
);

_sm_sqrt_gstar = CubicSplineInterpolation(
    LinRange(-4.7, 4.1, 276),
    _sm_data.SqrtGstar,
    extrapolation_bc = Flat(),
);

sm_heff(T::Real) = _sm_heff(log10(T))
sm_geff(T::Real) = _sm_geff(log10(T))
sm_sqrt_gstar(T::Real) = _sm_sqrt_gstar(log10(T))

sm_entropy_density(T::Real) = 2π^2 / 45 * _sm_heff(log10(T)) * T^3
sm_energy_density(T::Real) = π^2 / 30 * _sm_geff(log10(T)) * T^4
