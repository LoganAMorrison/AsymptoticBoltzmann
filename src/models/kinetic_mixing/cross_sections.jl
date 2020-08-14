function cross_section_xx_νν(Q::Real, model::KineticMixing)
    Q <= 2model.mx && return zero(typeof(Q))
    (
        ALPHA_EM *
        model.eps^2 *
        model.gvxx^2 *
        sqrt(Q^2) *
        (2 * model.mx^2 + Q^2)
    ) / (
        24 *
        COS_THETA_WEAK^2 *
        sqrt(-4 * model.mx^2 + Q^2) *
        ((model.mv^2 - Q^2)^2 + model.mv^2 * model.width_v^2)
    )
end

function cross_section_xx_ll(Q::Real, ml::Float64, model::KineticMixing)
    (Q <= 2model.mx || Q <= 2ml) && return zero(typeof(Q))
    (
        ALPHA_EM *
        model.eps^2 *
        model.gvxx^2 *
        (2 * model.mx^2 + Q^2) *
        sqrt(-4 * ml^2 + Q^2) *
        (7 * ml^2 + 5 * Q^2)
    ) / (
        24 *
        COS_THETA_WEAK^2 *
        Q^2 *
        sqrt(-4 * model.mx^2 + Q^2) *
        ((model.mv^2 - Q^2)^2 + model.mv^2 * model.width_v^2)
    )
end

cross_section_xx_ee(Q::Real, model::KineticMixing) =
    cross_section_xx_ll(Q, ELECTRON_MASS, model)
cross_section_xx_μμ(Q::Real, model::KineticMixing) =
    cross_section_xx_ll(Q, MUON_MASS, model)
cross_section_xx_ττ(Q::Real, model::KineticMixing) =
    cross_section_xx_ll(Q, TAU_MASS, model)

function cross_section_xx_ququ(Q::Real, mu::Float64, model::KineticMixing)
    (Q <= 2model.mx || Q <= 2mu) && return zero(typeof(Q))
    (
        ALPHA_EM *
        model.eps^2 *
        model.gvxx^2 *
        (2 * model.mx^2 + Q^2) *
        sqrt(-4 * mu^2 + Q^2) *
        (7 * mu^2 + 17 * Q^2)
    ) / (
        72 *
        COS_THETA_WEAK^2 *
        Q^2 *
        sqrt(-4 * model.mx^2 + Q^2) *
        ((model.mv^2 - Q^2)^2 + model.mv^2 * model.width_v^2)
    )
end

cross_section_xx_uu(Q::Real, model::KineticMixing) =
    cross_section_xx_ququ(Q, UP_QUARK_MASS, model)
cross_section_xx_cc(Q::Real, model::KineticMixing) =
    cross_section_xx_ququ(Q, CHARM_QUARK_MASS, model)
cross_section_xx_tt(Q::Real, model::KineticMixing) =
    cross_section_xx_ququ(Q, TOP_QUARK_MASS, model)

function cross_section_xx_qdqd(Q::Real, md::Float64, model::KineticMixing)
    (Q <= 2model.mx || Q <= 2md) && return zero(typeof(Q))
    (
        ALPHA_EM *
        model.eps^2 *
        model.gvxx^2 *
        sqrt(-4 * md^2 + Q^2) *
        (2 * model.mx^2 + Q^2) *
        (-17 * md^2 + 5 * Q^2)
    ) / (
        72 *
        COS_THETA_WEAK^2 *
        Q^2 *
        sqrt(-4 * model.mx^2 + Q^2) *
        ((model.mv^2 - Q^2)^2 + model.mv^2 * model.width_v^2)
    )
end

cross_section_xx_dd(Q::Real, model::KineticMixing) =
    cross_section_xx_qdqd(Q, DOWN_QUARK_MASS, model)
cross_section_xx_ss(Q::Real, model::KineticMixing) =
    cross_section_xx_qdqd(Q, STRANGE_QUARK_MASS, model)
cross_section_xx_bb(Q::Real, model::KineticMixing) =
    cross_section_xx_qdqd(Q, BOTTOM_QUARK_MASS, model)

function cross_section_xx_hz(Q::Real, model::KineticMixing)
    (Q <= 2model.mx || Q <= HIGGS_MASS + Z_BOSON_MASS) && return zero(typeof(Q))
    (
        ALPHA_EM^2 *
        model.eps^2 *
        model.gvxx^2 *
        pi *
        sqrt(
            (
                (HIGGS_MASS - Z_BOSON_MASS - Q) *
                (HIGGS_MASS + Z_BOSON_MASS - Q) *
                (HIGGS_MASS - Z_BOSON_MASS + Q) *
                (HIGGS_MASS + Z_BOSON_MASS + Q)
            ) / (-4 * model.mx^2 + Q^2),
        ) *
        (2 * model.mx^2 + Q^2) *
        (
            HIGGS_MASS^4 + Z_BOSON_MASS^4 + 10 * Z_BOSON_MASS^2 * Q^2 + Q^4 -
            2 * HIGGS_MASS^2 * (Z_BOSON_MASS^2 + Q^2)
        ) *
        HIGGS_VEV^2
    ) / (
        48 *
        COS_THETA_WEAK^4 *
        Z_BOSON_MASS^2 *
        Q^5 *
        SIN_THETA_WEAK^2 *
        ((model.mv^2 - Q^2)^2 + model.mv^2 * model.width_v^2)
    )
end

function cross_section_xx_vv(Q::Real, model::KineticMixing)
    (Q <= 2model.mx || Q <= 2model.mv) && return zero(typeof(Q))
    (
        model.gvxx^4 * (
            (
                -2 *
                sqrt(-4 * model.mx^2 + Q^2) *
                sqrt(-4 * model.mv^2 + Q^2) *
                (4 * model.mx^4 + 2 * model.mv^4 + model.mx^2 * Q^2)
            ) / (model.mv^4 + model.mx^2 * (-4 * model.mv^2 + Q^2)) +
            (2 * model.mx^2 - 2 * model.mv^2 + Q^2) * log(
                (
                    -2 * model.mv^2 +
                    Q^2 +
                    sqrt(-4 * model.mx^2 + Q^2) * sqrt(-4 * model.mv^2 + Q^2)
                )^2 /
                (
                    2 * model.mv^2 - Q^2 +
                    sqrt(-4 * model.mx^2 + Q^2) * sqrt(-4 * model.mv^2 + Q^2)
                )^2,
            ) +
            (
                2 *
                (
                    -4 * model.mx^4 +
                    2 * model.mv^2 * Q^2 +
                    model.mx^2 * (-2 * model.mv^2 + Q^2)
                ) *
                log(
                    (
                        -2 * model.mv^2 +
                        Q^2 +
                        sqrt(-4 * model.mx^2 + Q^2) *
                        sqrt(-4 * model.mv^2 + Q^2)
                    )^2 /
                    (
                        2 * model.mv^2 - Q^2 +
                        sqrt(-4 * model.mx^2 + Q^2) *
                        sqrt(-4 * model.mv^2 + Q^2)
                    )^2,
                )
            ) / (-2 * model.mv^2 + Q^2)
        )
    ) / (16 * pi * Q^2 * (-4 * model.mx^2 + Q^2))
end

function cross_section_xx_tot(Q::Real, model::KineticMixing)
    (
        3cross_section_xx_νν(Q, model) +
        cross_section_xx_ee(Q, model) +
        cross_section_xx_μμ(Q, model) +
        cross_section_xx_ττ(Q, model) +
        cross_section_xx_uu(Q, model) +
        cross_section_xx_cc(Q, model) +
        cross_section_xx_tt(Q, model) +
        cross_section_xx_dd(Q, model) +
        cross_section_xx_ss(Q, model) +
        cross_section_xx_bb(Q, model) +
        cross_section_xx_hz(Q, model) +
        cross_section_xx_vv(Q, model)
    )
end


function thermal_cross_section_xx_tot(x::Float64, model::KineticMixing)
    den = 2 * besselkx(2, x)
    pf = x / den^2 / 2

    integrand(z::Real) = (
        cross_section_xx_tot(model.mx * z, model) *
        z^2 *
        (z^2 - 4) *
        besselkx(1, x * z) *
        exp(-x * (z - 2))
    )

    res_loc = model.mv / model.mx
    thr_loc = 2res_loc

    if thr_loc > 2
        # We will pass a threshold in this case
        if thr_loc > 2
            # We will also pass a resonance
            return pf * quadgk(integrand, 2, res_loc, thr_loc, Inf)[1]
        else
            return pf * quadgk(integrand, 2, thr_loc, Inf)[1]
        end
    else
        # We can go directly to vectors. No threshold
        return pf * quadgk(integrand, 2, Inf)[1]
    end
end

function thermal_cross_section_xx_tot(x::Real, model::KineticMixing)
    den = 2 * besselk(2, x)
    pf = x / den^2 / 2

    integrand(z::Real) = (
        cross_section_xx_tot(model.mx * z, model) *
        z^2 *
        (z^2 - 4) *
        besselk(1, x * z)
    )

    res_loc = model.mv / model.mx
    thr_loc = 2res_loc

    if thr_loc > 2
        # We will pass a threshold in this case
        if thr_loc > 2
            # We will also pass a resonance
            return pf * quadgk(integrand, 2, res_loc, thr_loc, Inf)[1]
        else
            return pf * quadgk(integrand, 2, thr_loc, Inf)[1]
        end
    else
        # We can go directly to vectors. No threshold
        return pf * quadgk(integrand, 2, Inf)[1]
    end
end
