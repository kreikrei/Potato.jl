module Potato

using JuMP
using XLSX
using DataFrames
using Distances
using Coluna
using GLPK

include("base.jl")
include("core.jl")

#base
export extract!
export V,dist,K,T,d

end
