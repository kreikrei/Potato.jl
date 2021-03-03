module Potato

using JuMP
using XLSX
using DataFrames
using Distances
using UUIDs
using Statistics

include("struct.jl")
include("settings.jl")
include("base.jl")
include("core.jl")

#struct
export vtx,veh,col,node,dv,β,bound,S

#settings
export set_optimizer!
export get_optimizer
export reset_optimizer
export set_slack_coeff!
export set_surp_coeff!
export su_C
export sl_C
export silence!
export silent

#base
export extract!
export V,dist,K,T,d
export root

#core
export passes
export callSub,buildSub!
export master,getDuals
export sub,getCols
export colGen
export origin
export Q,s,f
export separate

end
