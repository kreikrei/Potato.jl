using Revise
using JuMP
using Gurobi
using Potato

path = joinpath(@__DIR__,"abs1n20.xlsx")
extract!(path;f="euclidean")

path = joinpath(DEPOT_PATH[1],"dev","Potato\\test\\A1.xlsx")
extract!(path;f="haversine")

GUROBI_ENV = Gurobi.Env()
set_optimizer!(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
set_slack_coeff!(1000.0)
set_surp_coeff!(-1000.0)

prioritize_vertex!()
prioritize_vehicle!()

test_root = root()

#=test_master = master(test_root)
println(test_master)
test_duals = getDuals(test_master)
test_sub = sub(test_root,test_duals)
test_cols = getCols(test_sub)

push!(test_root.columns,test_cols)=#
@time colGen(test_root,maxCG=Inf,track=true)

origin(test_root)

R = Dict(1:length(test_root.columns) .=> test_root.columns)
θ = value.(master(test_root).obj_dict[:θ])
new = find_branch(R,θ)
println(new)
s(new,R,θ)
to_push = bound(new,"<=",0)
push!(test_root.bounds,to_push)

test_root = root()
@time res = leaf(test_root,Inf,1000.0)

res.upperBound
tes = filter(p -> p.status[end] == "NEW_BOUND",res.visited)

test_master = master(tes[end])
origin(tes[end])

traverse(1000.0)

tes_I = value.(test_master.obj_dict[:I]).data
h_i = [V(i).h for i in V()]

sum(tes_I .* h_i)
dist(1,5) + dist(5,4) + dist(4,1)



[V(i).START for i in V()]

[d(i) for i in V()]


coluna = optimizer_with_attributes(
    Coluna.Optimizer,
    "params" => Coluna.Params(solver = Coluna.Algorithm.TreeSearchAlgorithm()),
    "default_optimizer" => Gurobi.Optimizer
)

colindex = Vector{NamedTuple}()
for k in K(), t in T(), m in 1:sum(K(k).BP[i] for i in K(k).cover)
    push!(colindex,(k=k,t=t,m=m))
end

@axis(𝕂𝕋𝕄, colindex)



𝕂𝕋𝕄[2].indice.k

𝕄 = Dict()
for k in K()
    𝕄[k] = collect(1:sum(K(k).BP[i] for i in K(k).cover))
end


@axis(𝕋, T())

ori = BlockModel(coluna)
@variable(ori, I[V(), vcat(first(T()) - 1, T())])

@variable(ori, u[k = 𝕂, t = 𝕋, m = 𝕄[k], i = V()] >= 0, Int)
@variable(ori, v[k = 𝕂, t = 𝕋, m = 𝕄[k], i = V()] >= 0, Int)
@variable(ori, l[k = 𝕂, t = 𝕋, m = 𝕄[k], i = V(), j = V()] >= 0, Int)
@variable(ori, o[k = 𝕂, t = 𝕋, m = 𝕄[k], i = V()], Bin)
@variable(ori, x[k = 𝕂, t = 𝕋, m = 𝕄[k], i = V(), j = V()], Bin)

#BASIC CONSTRAINTS
@constraint(ori, [k = 𝕂, t = 𝕋, m = 𝕄[k]],
    sum(u[k,t,m,i] for i in K(k).cover) - sum(v[k,t,m,i] for i in K(k).cover) == 0
)

@constraint(ori, [k = 𝕂, t = 𝕋, m = 𝕄[k], i = K(k).cover],
    sum(l[k,t,m,j,i] for j in K(k).cover) - sum(l[k,t,m,i,j] for j in K(k).cover) ==
    u[k,t,m,i] - v[k,t,m,i]
)

@constraint(ori, [k = 𝕂, t = 𝕋, m = 𝕄[k], i = K(k).cover],
    sum(x[k,t,m,j,i] for j in K(k).cover) - sum(x[k,t,m,i,j] for j in K(k).cover) == 0
)

@constraint(ori, [k = 𝕂, t = 𝕋, m = 𝕄[k], i = K(k).cover],
    v[k,t,m,i] <= K(k).Q * o[k,t,m,i]
) #VZ

@constraint(ori, [k = 𝕂, t = 𝕋, m = 𝕄[k], i = K(k).cover, j = K(k).cover],
    l[k,t,m,i,j] <= K(k).Q * x[k,t,m,i,j]
) #XL

@constraint(ori, [k = 𝕂,t = 𝕋, m = 𝕄[k]],
    sum(o[k,t,m,i] for i in K(k).cover) <= 1
) #one start

@constraint(ori, [k = 𝕂, t = 𝕋, i = K(k).cover],
    sum(o[k,t,m,i] for m in 𝕄[k]) <= K(k).BP[i]
)

@objective(ori, Min,
    sum(
        sum(
            dist(i,j) * (
                K(k).vx * x[k,t,m,i,j] +
                K(k).vl * l[k,t,m,i,j]
            )
            for i in K(k).cover, j in K(k).cover
        ) +
        sum(
            K(k).fd * u[k,t,m,i]
            for i in K(k).cover
        ) +
        sum(
            K(k).fp * o[k,t,m,i]
            for i in K(k).cover
        )
        for k in 𝕂, t in 𝕋, m in 𝕄[k]
    ) + #column costs
    sum(
        V(i).h * I[i,t]
        for i in V(), t in T()
    )
)

@constraint(ori, [i = V(), t = T()],
    I[i,t - 1] +
    sum(u[k,t,m,i] for k in 𝕂, m in 𝕄[k]) ==
    sum(v[k,t,m,i] for k in 𝕂, m in 𝕄[k]) +
    d(i,t) + I[i,t]
)

@constraint(ori, [i = V(), t = T()],
    V(i).MIN <= I[i,t] <= V(i).MAX #inventory capacity interval
)

@constraint(ori, [i = V()],
    I[i,first(T())-1] == V(i).START #starting inventory level
)

@dantzig_wolfe_decomposition(ori, decomposition, 𝕂)

master = getmaster(decomposition)
subproblems = getsubproblems(decomposition)

for k in K()
    specify!(
        subproblems[k], lower_multiplicity = 0,
        upper_multiplicity = sum(M[k]) * length(T())
    )
end

optimize!(ori)

value.(u)

println(value.(u))
println("+++++++++++++++")
