abstract type AbstractDMModel end

dm_mass(::AbstractDMModel) = throw("'dm_mass' isn't implemented for this model")
dm_dof(::AbstractDMModel) = throw("'dm_dof' isn't implemented for this model")
dm_spin2(::AbstractDMModel) =
    throw("'dm_spin2' isn't implemented for this model")
dm_annihilation_cross_section(Q::Real, ::AbstractDMModel) =
    throw("'dm_annihilation_cross_section' isn't implemented for this model")
dm_thermal_cross_section(x::Real, ::AbstractDMModel) =
    throw("'dm_thermal_cross_section' isn't implemented for this model")
