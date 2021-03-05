# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

passes(i) = [k for k in K() if (i in K(k).cover)]
s(key::S,R::Dict,θ) = sum(θ[q.r,q.k,q.t] for q in Q(key,R))
issinteger(val) = abs(round(val) - val) < 1e-8

function Q(key::S,R::Dict) #key pasti adalah S
    if isempty(key.seq)
        #if formnya S(k,β[])
        q = Vector{NamedTuple}()
        for r in keys(R)
            push!(q,(r=r,k=key.k,t=key.t))
        end
        return q
    else
        #if formnya S(k,β[β(i,v),...])
        q = Vector{Vector{NamedTuple}}()
        for tes in key.seq
            res = Vector{NamedTuple}()
            for r in keys(R)
                if R[r][key.k,key.t].u[tes.i] >= tes.v
                    push!(res,(r=r,k=key.k,t=key.t))
                end
            end
            push!(q,res) #add col each bound
        end
        return reduce(intersect,q)
    end
    #column class extraction
end

const vehicle_priority = Ref{Any}(nothing)
k_priority() = vehicle_priority[] #clear

const vertex_priority = Ref{Any}(nothing)
i_priority() = vertex_priority[] #clear

function prioritize_vertex!()
    #determine branching priorities
    res = JuMP.Containers.DenseAxisArray{Float64}(undef,V())
    res .= 0
    for i in V()
        slack = (V(i).START - sum(d(i,t) for t in T())) - V(i).MIN #after demand
        surplus = V(i).MAX - (V(i).START - sum(d(i,t) for t in T())) #after demand
        spare = min(slack,surplus)/(V(i).MAX - V(i).MIN)
        res[i] = 1/abs(spare)
    end
    return vertex_priority[] = res
end #clear

function i_priority(cover)
    if isa(cover,Int64)
        return vertex_priority[][cover] #bentuknya DenseAxisArray
    else
        new = Dict()
        for i in cover
            new[i] = i_priority(i)
        end
        return new
    end
end

function prioritize_vehicle!()
    cov = DataFrame(k=Int64[],cov=Int64[])
    for k in K()
        append!(cov,DataFrame(k=k,cov=length(K(k).cover)))
    end
    sort!(cov,:cov,rev=true)

    return vehicle_priority[] = cov.k
end

function find_branch(R,θ)
    #iterate over k_priority()
    for k in k_priority(), t in T()
        key = S(k,t,β[])
        if !issinteger(s(key,R,θ))
            return key
        end
    end

    for k in k_priority(),t in T()
        key = S(k,t,β[])
        rec = separate(key,deepcopy(K(k).cover),[],R,θ)
        if !isempty(rec) #rec is a vector of S
            lastcomp = [candidate.seq[end].i for candidate in rec]
            to_collect = findmax(i_priority(lastcomp))[2]
            to_return = filter(p -> p.seq[end].i == to_collect,rec)
            return to_return[1]
        end
    end
end

function separate(Seq::S,cover,record,R,θ)
    F = fractional_column(Seq,R,θ)

    if isempty(F)
        return record
    end

    found = false
    for i in cover
        v = find_separator(i,F,R)
        newSeq = S(Seq.k,Seq.t,β[β(i,v)])
        if !issinteger(s(newSeq,R,θ)) > 0
            push!(record,newSeq)
            found = true
        end
    end
    if found
        return record
    end

    i✶ = findmax(i_priority(cover))[2] #determine next highest priority

    v = find_separator(i✶,F,R)
    newSeq = S(Seq.k,Seq.t,β[β(i✶,v)]) #BUILD NEWSEQ

    to_remove = findall(p -> p == i✶,cover)
    splice!(cover,to_remove) #BUILDNEWCOVER

    record = separate(newSeq,cover,record,R,θ)
    return record
end

function fractional_column(Seq::S,R::Dict,θ)
    F = Vector{NamedTuple}()
    for r in θ.axes[1], k in θ.axes[2], t in θ.axes[3]
        if !issinteger(θ[r,k,t])
            push!(F,(r=r,k=k,t=t))
        end
    end

    return intersect(F,Q(Seq,R))
end

function find_separator(i::Int64,F,R)
    test = Vector{Int64}()
    for f in F
        push!(test,round(R[f.r][(f.k,f.t)].u[i]))
    end
    unique!(test)

    return ceil(median(test))
end

function f(key::S,R::Dict,θ)
    return sum(θ[q.r,q.k,q.t] - floor(θ[q.r,q.k,q.t]) for q in Q(key,R))
end

const subproblems = Ref{Any}(nothing)
callSub() = subproblems[]
callSub(k,t) = subproblems[][(k,t)]

function buildSub!(n::node)
    R = Dict{Tuple,Model}()

    @sync begin
        @inbounds for k in K(), t in T()
            @async begin
                sp = Model(get_optimizer())
                set_silent(sp)

                if solver_name(sp) == "Gurobi"
                    set_optimizer_attribute(sp,"MIPFocus",2)
                    set_optimizer_attribute(sp,"NodefileStart",0.5)
                    set_optimizer_attribute(sp, "NumericFocus",3)
                end

                #VARIABLE DEFINITION
                @variable(sp, u[K(k).cover] >= 0, Int) #peti diantar
                @variable(sp, v[K(k).cover] >= 0, Int) #peti dijemput
                @variable(sp, l[K(k).cover, K(k).cover] >= 0, Int) #load kendaraan
                @variable(sp, o[i = K(k).cover] <= K(k).BP[i], Bin) #origin
                @variable(sp, x[K(k).cover, K(k).cover], Bin) #penggunaan segmen

                @constraint(sp,
                    sum(u[i] for i in K(k).cover) - sum(v[i] for i in K(k).cover) == 0
                ) #all pickup delivered
                @constraint(sp, [i = K(k).cover],
                    sum(l[j,i] for j in K(k).cover) - sum(l[i,j] for j in K(k).cover) ==
                    u[i] - v[i]
                ) #load balance
                @constraint(sp,
                    [i = K(k).cover, j = K(k).cover], l[i,j] <= K(k).Q * x[i,j]
                ) #XL
                @constraint(sp, [i = K(k).cover],
                    sum(x[j,i] for j in K(k).cover) - sum(x[i,j] for j in K(k).cover) == 0
                ) #traversal
                @constraint(sp, [i = K(k).cover], v[i] <= K(k).Q * o[i]) #VZ
                @constraint(sp, sum(o[i] for i in K(k).cover) <= 1) #one start

                F = Dict(1:length(n.bounds) .=> n.bounds)
                uB = filter(f -> last(f).sense=="<="&&last(f).S.k==k && last(f).S.t == t, F)
                lB = filter(f -> last(f).sense==">="&&last(f).S.k==k && last(f).S.t == t, F)

                @variable(sp, g[keys(uB)], Bin)
                @variable(sp, h[keys(lB)], Bin)

                q = col(u,v,l,o,x)

                for j in keys(uB)
                    η = @variable(sp, [F[j].S.seq], Bin)
                    @constraint(sp, g[j] >= 1 - sum((1 - η[e]) for e in F[j].S.seq))
                    @constraint(sp, [e = F[j].S.seq],
                    (K(k).Q - e.v + 1) * η[e] >= getproperty(q,:u)[e.i] - e.v + 1)
                end

                for j in keys(lB)
                    η = @variable(sp, [F[j].S.seq], Bin)
                    @constraint(sp, [e = F[j].S.seq], h[j] <= η[e])
                    @constraint(sp, [e = F[j].S.seq],
                    e.v * η[e] <= getproperty(q,:u)[e.i])
                end

                optimize!(sp)

                R[(k,t)] = sp
            end
        end
    end

    return subproblems[] = R
end

function master(n::node)
    mp = Model(get_optimizer())
    set_silent(mp)

    if solver_name(mp) == "Gurobi"
        set_optimizer_attribute(mp,"MIPFocus",2)
        set_optimizer_attribute(mp,"NodefileStart",0.5)
        set_optimizer_attribute(mp, "NumericFocus",3)
    end

    R = Dict(1:length(n.columns) .=> n.columns)

    # ================================
    #    MODEL CONSTRUCTION
    # ================================
    @variable(mp, θ[keys(R), K(), T()] >= 0)
    @variable(mp, I[i = V(), vcat(first(T()) - 1, T())])
    @variable(mp, 0 <= slack[k = K(), t = T(), i = K(k).cover] <= n.stab.slLim[(k,t)])
    @variable(mp, 0 <= surp[k = K(), t = T(), i = K(k).cover] <= n.stab.suLim[(k,t)])

    @objective(mp, Min,
        sum(
            θ[r,k,t] * (
                sum(
                    dist(i,j) * (
                        K(k).vx * R[r][(k,t)].x[i,j] +
                        K(k).vl * R[r][(k,t)].l[i,j]
                    )
                    for i in K(k).cover, j in K(k).cover
                ) +
                sum(
                    K(k).fd * R[r][(k,t)].u[i]
                    for i in K(k).cover
                ) +
                sum(
                    K(k).fp * R[r][(k,t)].o[i]
                    for i in K(k).cover
                )
            )
            for r in keys(R), k in K(), t in T()
        ) + #column costs
        sum(
            V(i).h * I[i,t]
            for i in V(), t in T()
        ) + #inventory costs
        sum(
            n.stab.slCoeff * slack[k,t,i]
            for k in K(), t in T(), i in K(k).cover
        ) - #stabilizer
        sum(
            n.stab.suCoeff * surp[k,t,i]
            for k in K(), t in T(), i in K(k).cover
        ) #stabilizer
    )

    @constraint(mp, λ[i = V(), t = T()],
        I[i,t-1] + sum(R[r][(k,t)].u[i] * θ[r,k,t] for r in keys(R), k in passes(i)) +
        sum(slack[k,t,i] for k in passes(i))  ==
        sum(R[r][(k,t)].v[i] * θ[r,k,t] for r in keys(R), k in passes(i)) +
        sum(surp[k,t,i] for k in passes(i)) + d(i,t) + I[i,t]
    ) #inventory balance

    @constraint(mp, δ[k = K(), t = T(), i = K(k).cover],
        sum(R[r][(k,t)].o[i] * θ[r,k,t] for r in keys(R)) <= K(k).BP[i]
    ) #limit starting point

    @constraint(mp, ϵ[k = K(), t = T()],
        sum(θ[r,k,t] for r in keys(R)) <= sum(K(k).BP[i] for i in K(k).cover)
    ) #multiplicity

    @constraint(mp, [i = V(), t = T()],
        V(i).MIN <= I[i,t] <= V(i).MAX
    ) #inventory max min

    @constraint(mp, [i = V()],
        I[i,first(T())-1] == V(i).START
    ) #inventory start

    F = Dict(1:length(n.bounds) .=> n.bounds)
    uB = filter(f -> last(f).sense == "<=",F)
    lB = filter(f -> last(f).sense == ">=",F)

    @constraint(mp, ρ[j = keys(uB)], s(F[j].S,R,θ) <= F[j].κ)
    @constraint(mp, σ[j = keys(lB)], s(F[j].S,R,θ) >= F[j].κ)

    optimize!(mp)

    return mp
end

function getDuals(mp::Model)
    λ = dual.(mp.obj_dict[:λ])
    δ = dual.(mp.obj_dict[:δ])
    ϵ = dual.(mp.obj_dict[:ϵ])

    ρ = dual.(mp.obj_dict[:ρ])
    σ = dual.(mp.obj_dict[:σ])

    return dv(λ,δ,ϵ,ρ,σ)
end

function sub(n::node,duals::dv)
    @sync begin
        @inbounds for k in K(), t in T()
            @async begin
                sp = callSub(k,t)

                F = Dict(1:length(n.bounds) .=> n.bounds)
                uB = filter(f->last(f).sense == "<="&&last(f).S.k == k&&last(f).S.t == t, F)
                lB = filter(f->last(f).sense == ">="&&last(f).S.k == k&&last(f).S.t == t, F)

                #ADD OBJECTIVE
                @objective(sp, Min,
                    sum(
                        dist(i,j) * (
                            K(k).vx * sp.obj_dict[:x][i,j] +
                            K(k).vl * sp.obj_dict[:l][i,j]
                        )
                        for i in K(k).cover, j in K(k).cover
                    ) +
                    sum(K(k).fd * sp.obj_dict[:u][i] for i in K(k).cover) +
                    sum(K(k).fp * sp.obj_dict[:o][i] for i in K(k).cover) -
                    sum(
                        (sp.obj_dict[:u][i] - sp.obj_dict[:v][i]) * duals.λ[i,t]
                        for i in K(k).cover
                    ) -
                    sum(sp.obj_dict[:o][i] * duals.δ[k,t,i] for i in K(k).cover) -
                    duals.ϵ[k,t] -
                    sum(sp.obj_dict[:g][j] * duals.ρ[j] for j in keys(uB)) -
                    sum(sp.obj_dict[:h][j] * duals.σ[j] for j in keys(lB))
                )

                optimize!(sp)
                #println("($k): $(objective_value(sp))")
            end
        end
    end

    return callSub()
end

function getCols(sp)
    if isa(sp,Model)
        u = value.(sp.obj_dict[:u])
        v = value.(sp.obj_dict[:v])
        l = value.(sp.obj_dict[:l])
        o = value.(sp.obj_dict[:o])
        x = value.(sp.obj_dict[:x])

        return col(u,v,l,o,x)
    elseif isa(sp,Dict)
        new = Dict{Tuple,col}()
        for r in keys(sp)
            new[r] = getCols(sp[r])
        end

        return new
    end
end

function colvals()
    collection = 0
    for k in K(), t in T()
        collection += objective_value(callSub(k,t))
    end

    return sum(collection)
end

function updateStab!(stab::stabilizer,param::Float64)
    for kt in keys(stab.slLim)
        stab.slLim[kt] = floor(param * stab.slLim[kt])
    end
    for kt in keys(stab.suLim)
        stab.suLim[kt] = floor(param * stab.suLim[kt])
    end

    return stab
end

checkStab(mp::Model) = sum(value.(mp.obj_dict[:slack])) + sum(value.(mp.obj_dict[:surp]))

function colGen(n::node;maxCG::Float64,track::Bool)
    terminate = false
    iter = 0
    buildSub!(n)

    while !terminate
        if iter < maxCG
            mp = master(n)

            if has_values(mp) && has_duals(mp)
                if track #print master problem obj
                    println("obj: $(objective_value(mp))")
                end

                duals = getDuals(mp)
                sp = sub(n,duals)

                if track #print subproblem price
                    println("price: $(colvals())")
                end

                if isapprox(colvals(),0,atol = 1e-8) || colvals() > 0
                    if isapprox(checkStab(mp),0,atol = 1e-8)
                        terminate = true #action
                        push!(n.status,"EVALUATED") #report
                        if track
                            println("EVALUATED")
                        end
                    else
                        updateStab!(n.stab,0.2) #action
                        push!(n.status,"STABILIZED") #report
                        if track
                            println("STABILIZED")
                        end
                    end
                else
                    push!(n.columns,getCols(sp)) #action
                    push!(n.status,"ADD_COLUMN") #report
                    if track
                        println("ADD_COLUMN")
                    end
                end

                iter += 1 #iteration update
            else
                terminate = true #action
                push!(n.status,"NO_SOLUTION")
                if track
                    println("NO_SOLUTION")
                end
            end
        else
            terminate = true #action
            pop!(n.columns)
            updateStab!(n.stab,1e-8)
            push!(n.status,"EVALUATED") #report
            if track
                println("EVALUATED")
            end
        end
    end

    if n.status[end] == "NO_SOLUTION"
        println("NODE $(n.self) FAILED.")
    else
        if integerCheck(n)
            push!(n.status,"INTEGER")
            println("NODE $(n.self) INTEGER")
        else
            println("NODE $(n.self) FINISHED.")
        end
    end

    return n
end

function origin(n::node)
    o = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),K(),T())
    u = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),K(),T())
    v = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),K(),T())
    x = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),V(),K(),T())
    l = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),V(),K(),T())

    o .= 0
    u .= 0
    v .= 0
    l .= 0
    x .= 0

    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    for k in K(), t in T(), i in K(k).cover
        o[i,k,t] = sum(R[r][(k,t)].o[i] * θ[r,k,t] for r in keys(R))
        u[i,k,t] = sum(R[r][(k,t)].u[i] * θ[r,k,t] for r in keys(R))
        v[i,k,t] = sum(R[r][(k,t)].v[i] * θ[r,k,t] for r in keys(R))
    end

    for k in K(), t in T(), i in K(k).cover, j in K(k).cover
        x[i,j,k,t] = sum(R[r][(k,t)].x[i,j] * θ[r,k,t] for r in keys(R))
        l[i,j,k,t] = sum(R[r][(k,t)].l[i,j] * θ[r,k,t] for r in keys(R))
    end

    return col(u,v,l,o,x)
end
