module Potato

using JuMP
using XLSX
using DataFrames
using Distances
using Coluna
using GLPK

include("base.jl")
include("core.jl")
include("struct.jl")

#struct
export vtx,veh

#base
export extract!
export V,dist,K,T,d

end
