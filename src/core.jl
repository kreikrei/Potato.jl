# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

passes(i) = [k for k in K() if (i in K(k).cover)]
s(key,R,θ) = sum(θ[q.r,q.k] for q in Q(key,R))
issinteger(val) = abs(round(val) - val) < 1e-8

function Q(key,R) #key pasti adalah S
    if isempty(key.seq)
        #if formnya S(k,β[])
        q = Vector{NamedTuple}()
        for r in keys(R)
            push!(q,(r=r,k=key.k))
        end
        return q
    else
        #if formnya S(k,β[β(i,v),...])
        q = Vector{Vector{NamedTuple}}()
        for tes in key.seq
            res = Vector{NamedTuple}()
            for r in keys(R)
                if R[r][key.k].u[tes.i] >= tes.v
                    push!(res,(r=r,k=key.k))
                end
            end
            push!(q,res) #add col each bound
        end
        return reduce(intersect,q)
    end
    #column class extraction
end

function separate(R,θ)
    #separation of fractional column
end

function f(key,R,θ)
    #find fractionality of key (S)
    return sum(θ[q.r,q.k] - floor(θ[q.r,q.k]) for q in Q(key,R))
end

const subproblems = Ref{Any}(nothing)
callSub() = subproblems[]
callSub(k) = subproblems[][k]

function buildSub!(n::node)
    R = Dict{Int64,Model}()

    @inbounds for k in K()
        sp = Model(get_optimizer())
        set_silent(sp)

        if solver_name(sp) == "Gurobi"
            set_optimizer_attribute(sp,"MIPFocus",2)
            set_optimizer_attribute(sp,"NodefileStart",0.5)
            set_optimizer_attribute(sp, "NumericFocus",3)
        end

        #VARIABLE DEFINITION
        @variable(sp, u[K(k).cover] >= 0, Int)
        @variable(sp, v[K(k).cover] >= 0, Int)
        @variable(sp, l[K(k).cover, K(k).cover] >= 0, Int)
        @variable(sp, o[i = K(k).cover] <= K(k).BP[i], Bin)
        @variable(sp, x[K(k).cover, K(k).cover], Bin)

        @constraint(sp,
            sum(u[i] for i in K(k).cover) - sum(v[i] for i in K(k).cover) == 0
        ) #all pickup delivered
        @constraint(sp, [i = K(k).cover],
            sum(l[j,i] for j in K(k).cover) - sum(l[i,j] for j in K(k).cover) ==
            u[i] - v[i]
        ) #load balance
        @constraint(sp, [i = K(k).cover, j = K(k).cover], l[i,j] <= K(k).Q * x[i,j]) #XL
        @constraint(sp, [i = K(k).cover],
            sum(x[j,i] for j in K(k).cover) - sum(x[i,j] for j in K(k).cover) == 0
        ) #traversal
        @constraint(sp, [i = K(k).cover], v[i] <= K(k).Q * o[i]) #VZ
        @constraint(sp, sum(o[i] for i in K(k).cover) <= 1) #one start

        F = Dict(1:length(n.bounds) .=> n.bounds)
        uB = filter(f -> last(f).sense == "<=" && last(f).S.k == k, F) #k of S main di sini
        lB = filter(f -> last(f).sense == ">=" && last(f).S.k == k, F) #k of S main di sini

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

        R[k] = sp
    end

    subproblems[] = R
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
    @variable(mp, θ[keys(R), K()] >= 0)
    @variable(mp, V(i).MIN <= I[i = V()] <= V(i).MAX)
    @variable(mp, 0 <= slack[k = K(), i = K(k).cover] <= n.stab.slLim[k])
    @variable(mp, 0 <= surp[k = K(), i = K(k).cover] <= n.stab.suLim[k])

    @objective(mp, Min,
        sum(
            θ[r,k] * (
                sum(
                    dist(i,j) * (
                        K(k).vx * R[r][k].x[i,j] +
                        K(k).vl * R[r][k].l[i,j]
                    )
                    for i in K(k).cover, j in K(k).cover
                ) +
                sum(
                    K(k).fd * R[r][k].u[i]
                    for i in K(k).cover
                ) +
                sum(
                    K(k).fp * R[r][k].o[i]
                    for i in K(k).cover
                )
            )
            for r in keys(R), k in K()
        ) + #column costs
        sum(
            V(i).h * I[i]
            for i in V()
        ) + #inventory costs
        sum(
            n.stab.slCoeff * slack[k,i]
            for k in K(), i in K(k).cover
        ) - #stabilizer
        sum(
            n.stab.suCoeff * surp[k,i]
            for k in K(), i in K(k).cover
        ) #stabilizer
    )

    @constraint(mp, λ[i = V()],
        V(i).START +
        sum(R[r][k].u[i] * θ[r,k] for r in keys(R), k in passes(i)) +
        sum(slack[k,i] for k in passes(i))  ==
        sum(R[r][k].v[i] * θ[r,k] for r in keys(R), k in passes(i)) +
        sum(surp[k,i] for k in passes(i)) +
        d(i) + I[i]
    ) #inventory balance

    @constraint(mp, δ[k = K(), i = K(k).cover],
        sum(R[r][k].o[i] * θ[r,k] for r in keys(R)) <= K(k).BP[i]
    ) #limit starting point

    @constraint(mp, ϵ[k = K()],
        sum(θ[r,k] for r in keys(R)) <= sum(K(k).BP[i] for i in K(k).cover)
    ) #multiplicity

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
    @inbounds for k in K()
        sp = callSub()[k]

        F = Dict(1:length(n.bounds) .=> n.bounds)
        uB = filter(f -> last(f).sense == "<=" && last(f).S.k == k, F) #k of S main di sini
        lB = filter(f -> last(f).sense == ">=" && last(f).S.k == k, F) #k of S main di sini

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
                (sp.obj_dict[:u][i] - sp.obj_dict[:v][i]) * duals.λ[i]
                for i in K(k).cover
            ) -
            sum(sp.obj_dict[:o][i] * duals.δ[k,i] for i in K(k).cover) -
            duals.ϵ[k] -
            sum(sp.obj_dict[:g][j] * duals.ρ[j] for j in keys(uB)) -
            sum(sp.obj_dict[:h][j] * duals.σ[j] for j in keys(lB))
        )

        optimize!(sp)
        #println("($k): $(objective_value(sp))")
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
        new = Dict{Int64,col}()
        for r in keys(sp)
            new[r] = getCols(sp[r])
        end

        return new
    end
end

function colvals()
    collection = 0
    for k in K()
        collection += objective_value(callSub(k))
    end

    return sum(collection)
end

function updateStab!(stab::stabilizer,param::Float64)
    for k in keys(stab.slLim)
        stab.slLim[k] = floor(param * stab.slLim[k])
    end
    for k in keys(stab.suLim)
        stab.suLim[k] = floor(param * stab.suLim[k])
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
                        updateStab!(n.stab,0.5) #action
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
            push!(n.status,"EVALUATED") #report
            if track
                println("EVALUATED")
            end
        end
    end

    if n.status[end] == "NO_SOLUTION"
        println("NODE $(n.self) FAILED.")
    else
        println("NODE $(n.self) FINISHED.")
    end

    return n
end

function origin(n::node)
    o = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),K())
    u = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),K())
    v = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),K())
    x = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),V(),K())
    l = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),V(),K())

    o .= 0
    u .= 0
    v .= 0
    l .= 0
    x .= 0

    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    for k in K(), i in K(k).cover
        o[i,k] = sum(R[r][k].o[i] * θ[r,k] for r in keys(R))
        u[i,k] = sum(R[r][k].u[i] * θ[r,k] for r in keys(R))
        v[i,k] = sum(R[r][k].v[i] * θ[r,k] for r in keys(R))
    end

    for k in K(), i in K(k).cover, j in K(k).cover
        x[i,j,k] = sum(R[r][k].x[i,j] * θ[r,k] for r in keys(R))
        l[i,j,k] = sum(R[r][k].l[i,j] * θ[r,k] for r in keys(R))
    end

    return col(u,v,l,o,x)
end
