module Potato

using JuMP
using XLSX
using DataFrames
using Distances

include("struct.jl")
include("base.jl")
include("core.jl")

#struct
export vtx,veh,col,node,dv

#base
export extract!
export V,dist,K,T,d
export root

#core
export passes
export callSub,buildSub!
export master

end
