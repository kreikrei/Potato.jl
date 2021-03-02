# =========================================================================
#    BASIC STRUCTURES
# =========================================================================

struct vtx
    #IDENTIFIERS
    name::String
    type::String

    #VALS
    x::Float64 #xcoor
    y::Float64 #ycoor
    MAX::Float64 #max inventory level
    MIN::Float64 #min inventory level
    START::Float64 #starting inventory level

    #COSTS
    h::Float64 #inv cost per unit
end

struct veh
    #IDENTIFIERS
    name::String
    type::String

    #CHAR
    cover::Vector{Int64}
    BP::Dict{Int64,Int64}
    Q::Int64

    #COSTS
    vx::Float64
    vl::Float64
    fp::Float64
    fd::Float64
end

struct col
    #quantity
    u::JuMP.Containers.DenseAxisArray
    v::JuMP.Containers.DenseAxisArray
    l::JuMP.Containers.DenseAxisArray

    #decision
    o::JuMP.Containers.DenseAxisArray
    x::JuMP.Containers.DenseAxisArray
end
