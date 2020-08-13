struct KineticMixing
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
