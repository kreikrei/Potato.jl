# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

passes(i) = [k for k in K() if (i in K(k).cover)]

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

        @constraint(sp, [i = K(k).cover],
            sum(x[j,i] for j in K(k).cover) - sum(x[i,j] for j in K(k).cover) == 0
        ) #traversal

        @constraint(sp, [i = K(k).cover], v[i] <= K(k).Q * o[i]) #VZ
        @constraint(sp, [i = K(k).cover, j = K(k).cover], l[i,j] <= K(k).Q * x[i,j]) #XL

        @constraint(sp, sum(o[i] for i in K(k).cover) <= 1) #one start

        optimize!(sp)

        R[k] = sp
    end

    subproblems[] = R
end

function master(n::node)
    mp = Model(get_optimizer())
    set_silent(mp)

    if solver_name(mp) == "Gurobi"
        set_optimizer_attribute(sp,"MIPFocus",2)
        set_optimizer_attribute(sp,"NodefileStart",0.5)
        set_optimizer_attribute(sp, "NumericFocus",3)
    end

    R = Dict(1:length(n.columns) .=> n.columns)

    # ================================
    #    MODEL CONSTRUCTION
    # ================================
    @variable(mp, θ[keys(R), K()] >= 0)
    @variable(mp, V(i).MIN <= I[i = V()] <= V(i).MAX)
    @variable(mp, 0 <= slack[i = V()] <= n.stab.slLim[i])
    @variable(mp, 0 <= surp[i = V()] <= n.stab.suLim[i])

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
            n.stab.slCoeff * slack[i]
            for i in V()
        ) - #stabilizer
        sum(
            n.stab.suCoeff * surp[i]
            for i in V()
        ) #stabilizer
    )

    @constraint(mp, λ[i = V()],
        V(i).START[i] + sum(R[r][k].u[i] * θ[r,k] for r in keys(R), k in passes(i)) +
        slack[i] - surp[i] == sum(R[r][k].v[i] * θ[r,k] for r in keys(R), k in passes(i)) +
        d(i) + I[i]
    ) #inventory balance

    @constraint(mp, δ[k = K(), i = K(k).cover],
        sum(R[r][k].o[i] * θ[r,k] for r in keys(R)) <= K(k).BP[i]
    ) #limit starting point

    @constraint(mp, ϵ[k = K()],
        sum(θ[r,k] for r in keys(R)) <= sum(K(k).BP[i] for i in K(k).cover)
    ) #multiplicity

    optimize!(mp)

    return mp
end
