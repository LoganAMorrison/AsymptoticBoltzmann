struct KineticMixing <: AbstractDMModel
    mx::Float64
    mv::Float64
    gvxx::Float64
    eps::Float64
    width_v::Float64
end

function KineticMixing(mx::Float64, mv::Float64, gvxx::Float64, eps::Float64)
    model = KineticMixing(mx, mv, gvxx, eps, 0.0)
    width = width_v(model)
    KineticMixing(mx, mv, gvxx, eps, width)
end

dm_mass(model::KineticMixing) = model.mx
dm_dof(::KineticMixing) = 2.0
dm_spin2(::KineticMixing) = 2

# =================================
# Partial widths of vector mediator
# =================================

"""
    width_v_xx(model::KineticMixing)

Compute the partial width of the vector mediator into dark matter.
"""
function width_v_xx(model::KineticMixing)
    2model.mx >= model.mv && return 0.0
    (
        model.gvxx^2 *
        sqrt(-4 * model.mx^2 + model.mv^2) *
        (2 * model.mx^2 + model.mv^2)
    ) / (12 * model.mv^2 * pi)
end

"""
    width_v_νν(model::KineticMixing)

Compute the partial width of the vector mediator into neutrinos.
"""
function width_v_νν(model::KineticMixing)
    (ALPHA_EM * model.eps^2 * model.mv) / (24 * COS_THETA_WEAK^2)
end

"""
    width_v_ll(ml::Float64, model::KineticMixing)

Compute the partial width of the vector mediator into generic leptons.
"""
function width_v_ll(ml::Float64, model::KineticMixing)
    2ml >= model.mv && return 0.0
    (
        ALPHA_EM *
        model.eps^2 *
        sqrt(-4 * ml^2 + model.mv^2) *
        (7 * ml^2 + 5 * model.mv^2)
    ) / (24 * COS_THETA_WEAK^2 * model.mv^2)
end

"""
    width_v_ee(model::KineticMixing)

Compute the partial width of the vector mediator into electrons.
"""
width_v_ee(model::KineticMixing) = width_v_ll(ELECTRON_MASS, model)

"""
    width_v_μμ(model::KineticMixing)

Compute the partial width of the vector mediator into muons.
"""
width_v_μμ(model::KineticMixing) = width_v_ll(MUON_MASS, model)

"""
    width_v_ττ(model::KineticMixing)

Compute the partial width of the vector mediator into taus.
"""
width_v_ττ(model::KineticMixing) = width_v_ll(TAU_MASS, model)

"""
    width_v_ququ(model::KineticMixing)

Compute the partial width of the vector mediator into generic up-type quarks.
"""
function width_v_ququ(mu::Float64, model::KineticMixing)
    2mu >= model.mv && return 0.0
    (
        ALPHA_EM *
        model.eps^2 *
        sqrt(-4 * mu^2 + model.mv^2) *
        (7 * mu^2 + 17 * model.mv^2)
    ) / (72 * COS_THETA_WEAK^2 * model.mv^2)
end

"""
    width_v_uu(model::KineticMixing)

Compute the partial width of the vector mediator into generic up quarks.
"""
width_v_uu(model::KineticMixing) = width_v_ququ(UP_QUARK_MASS, model)

"""
    width_v_cc(model::KineticMixing)

Compute the partial width of the vector mediator into generic charms quarks.
"""
width_v_cc(model::KineticMixing) = width_v_ququ(CHARM_QUARK_MASS, model)

"""
    width_v_tt(model::KineticMixing)

Compute the partial width of the vector mediator into generic top quarks.
"""
width_v_tt(model::KineticMixing) = width_v_ququ(TOP_QUARK_MASS, model)

"""
    width_v_qdqd(model::KineticMixing)

Compute the partial width of the vector mediator into generic down-type quarks.
"""
function width_v_qdqd(md::Float64, model::KineticMixing)
    2md >= model.mv && return 0.0
    (
        ALPHA_EM *
        model.eps^2 *
        sqrt(-4 * md^2 + model.mv^2) *
        (-17 * md^2 + 5 * model.mv^2)
    ) / (72 * COS_THETA_WEAK^2 * model.mv^2)
end

"""
    width_v_dd(model::KineticMixing)

Compute the partial width of the vector mediator into generic down quarks.
"""
width_v_dd(model::KineticMixing) = width_v_qdqd(DOWN_QUARK_MASS, model)

"""
    width_v_ss(model::KineticMixing)

Compute the partial width of the vector mediator into generic strange quarks.
"""
width_v_ss(model::KineticMixing) = width_v_qdqd(STRANGE_QUARK_MASS, model)

"""
    width_v_bb(model::KineticMixing)

Compute the partial width of the vector mediator into generic bottom quarks.
"""
width_v_bb(model::KineticMixing) = width_v_qdqd(BOTTOM_QUARK_MASS, model)

"""
    width_v_hz(model::KineticMixing)

Compute the partial width of the vector mediator into generic Higgs + Z.
"""
function width_v_hz(model::KineticMixing)
    HIGGS_MASS + Z_BOSON_MASS >= model.mv && return 0.0

    (
        ALPHA_EM^2 *
        model.eps^2 *
        sqrt(
            -HIGGS_MASS^2 +
            (HIGGS_MASS^2 + model.mv^2 - Z_BOSON_MASS^2)^2 / (4 * model.mv^2),
        ) *
        (
            2 * model.mv^2 * Z_BOSON_MASS^2 +
            (-HIGGS_MASS^2 + model.mv^2 + Z_BOSON_MASS^2)^2 / 4.0
        ) *
        pi *
        HIGGS_VEV^2
    ) / (6 * COS_THETA_WEAK^4 * model.mv^4 * Z_BOSON_MASS^2 * SIN_THETA_WEAK^2)
end

"""
    width_v(model::KineticMixing)

Compute the total width of the vector mediator.
"""
function width_v(model::KineticMixing)
    (
        width_v_xx(model) +
        width_v_ee(model) +
        width_v_μμ(model) +
        width_v_ττ(model) +
        width_v_uu(model) +
        width_v_cc(model) +
        width_v_tt(model) +
        width_v_dd(model) +
        width_v_ss(model) +
        width_v_tt(model) +
        width_v_hz(model) +
        3width_v_νν(model)
    )
end

# =========================================================================
# ==== Cross sections =====================================================
# =========================================================================

"""
    cross_section_xx_νν(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into neutrinos at
center-of-mass energy `Q`.
"""
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

"""
    cross_section_xx_ll(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into generic leptons at
center-of-mass energy `Q`.
"""
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

"""
    cross_section_xx_ee(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into electons at
center-of-mass energy `Q`.
"""
cross_section_xx_ee(Q::Real, model::KineticMixing) =
    cross_section_xx_ll(Q, ELECTRON_MASS, model)

"""
    cross_section_xx_μμ(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into muon at
center-of-mass energy `Q`.
"""
cross_section_xx_μμ(Q::Real, model::KineticMixing) =
    cross_section_xx_ll(Q, MUON_MASS, model)

"""
    cross_section_xx_ττ(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into taus at
center-of-mass energy `Q`.
"""
cross_section_xx_ττ(Q::Real, model::KineticMixing) =
    cross_section_xx_ll(Q, TAU_MASS, model)

"""
    cross_section_xx_ququ(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into generic up-type quarks at
center-of-mass energy `Q`.
"""
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

"""
    cross_section_xx_uu(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into up quarks at
center-of-mass energy `Q`.
"""
cross_section_xx_uu(Q::Real, model::KineticMixing) =
    cross_section_xx_ququ(Q, UP_QUARK_MASS, model)

"""
    cross_section_xx_cc(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into charm quarks at
center-of-mass energy `Q`.
"""
cross_section_xx_cc(Q::Real, model::KineticMixing) =
    cross_section_xx_ququ(Q, CHARM_QUARK_MASS, model)

"""
    cross_section_xx_tt(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into top quarks at
center-of-mass energy `Q`.
"""
cross_section_xx_tt(Q::Real, model::KineticMixing) =
    cross_section_xx_ququ(Q, TOP_QUARK_MASS, model)

"""
    cross_section_xx_qdqd(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into generic down-type quarks at
center-of-mass energy `Q`.
"""
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

"""
    cross_section_xx_dd(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into down quarks at
center-of-mass energy `Q`.
"""
cross_section_xx_dd(Q::Real, model::KineticMixing) =
    cross_section_xx_qdqd(Q, DOWN_QUARK_MASS, model)

"""
    cross_section_xx_ss(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into strange quarks at
center-of-mass energy `Q`.
"""
cross_section_xx_ss(Q::Real, model::KineticMixing) =
    cross_section_xx_qdqd(Q, STRANGE_QUARK_MASS, model)

"""
    cross_section_xx_bb(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into bottom quarks at
center-of-mass energy `Q`.
"""
cross_section_xx_bb(Q::Real, model::KineticMixing) =
    cross_section_xx_qdqd(Q, BOTTOM_QUARK_MASS, model)

"""
    cross_section_xx_hz(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into Higgs + Z at
center-of-mass energy `Q`.
"""
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

"""
    cross_section_xx_vv(Q::Real, model::KineticMixing)

Compute cross section for DM annihilating into vector mediator at
center-of-mass energy `Q`.
"""
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

"""
    dm_annihilation_cross_section(Q::Real, model::KineticMixing)

Compute total annihilation cross section for DM center-of-mass energy `Q`.
"""
function dm_annihilation_cross_section(Q::Real, model::KineticMixing)
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

"""
    dm_thermal_cross_section(Q::Float64, model::KineticMixing)

Compute thermally-averaged annihilation cross section for DM at a scaled
temperatue `x=mx/T`
"""
function dm_thermal_cross_section(x::Float64, model::KineticMixing)
    den = 2 * besselkx(2, x)
    pf = x / den^2 / 2

    integrand(z::Real) = (
        dm_annihilation_cross_section(model.mx * z, model) *
        z^2 *
        (z^2 - 4) *
        besselkx(1, x * z) *
        exp(-x * (z - 2))
    )

    return pf * quadgk(integrand, 2, Inf)[1]
end

"""
    dm_thermal_cross_section(Q::Real, model::KineticMixing)

Compute thermally-averaged annihilation cross section for DM at a scaled
temperatue `x=mx/T`
"""
function dm_thermal_cross_section(x::Real, model::KineticMixing)
    den = 2 * besselk(2, x)
    pf = x / den^2 / 2

    integrand(z::Real) = (
        dm_annihilation_cross_section(model.mx * z, model) *
        z^2 *
        (z^2 - 4) *
        besselk(1, x * z)
    )
    return pf * quadgk(integrand, 2, Inf)[1]
end
