function width_v_xx(model::KineticMixing)
    2model.mx >= model.mv && return 0.0
    (
        model.gvxx^2 *
        sqrt(-4 * model.mx^2 + model.mv^2) *
        (2 * model.mx^2 + model.mv^2)
    ) / (12 * model.mv^2 * pi)
end

function width_v_νν(model::KineticMixing)
    (ALPHA_EM * model.eps^2 * model.mv) / (24 * COS_THETA_WEAK^2)
end

# Partial widths of VM to leptons

function width_v_ll(ml::Float64, model::KineticMixing)
    2ml >= model.mv && return 0.0
    (
        ALPHA_EM *
        model.eps^2 *
        sqrt(-4 * ml^2 + model.mv^2) *
        (7 * ml^2 + 5 * model.mv^2)
    ) / (24 * COS_THETA_WEAK^2 * model.mv^2)
end

width_v_ee(model::KineticMixing) = width_v_ll(ELECTRON_MASS, model)
width_v_μμ(model::KineticMixing) = width_v_ll(MUON_MASS, model)
width_v_ττ(model::KineticMixing) = width_v_ll(TAU_MASS, model)

# Partial widths of VM to up-type quarks

function width_v_ququ(mu::Float64, model::KineticMixing)
    2mu >= model.mv && return 0.0
    (
        ALPHA_EM *
        model.eps^2 *
        sqrt(-4 * mu^2 + model.mv^2) *
        (7 * mu^2 + 17 * model.mv^2)
    ) / (72 * COS_THETA_WEAK^2 * model.mv^2)
end

width_v_uu(model::KineticMixing) = width_v_ququ(UP_QUARK_MASS, model)
width_v_cc(model::KineticMixing) = width_v_ququ(CHARM_QUARK_MASS, model)
width_v_tt(model::KineticMixing) = width_v_ququ(TOP_QUARK_MASS, model)

# Partial widths of VM to down-type quarks

function width_v_qdqd(md::Float64, model::KineticMixing)
    2md >= model.mv && return 0.0
    (
        ALPHA_EM *
        model.eps^2 *
        sqrt(-4 * md^2 + model.mv^2) *
        (-17 * md^2 + 5 * model.mv^2)
    ) / (72 * COS_THETA_WEAK^2 * model.mv^2)
end

width_v_dd(model::KineticMixing) = width_v_qdqd(DOWN_QUARK_MASS, model)
width_v_ss(model::KineticMixing) = width_v_qdqd(STRANGE_QUARK_MASS, model)
width_v_bb(model::KineticMixing) = width_v_qdqd(BOTTOM_QUARK_MASS, model)

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
