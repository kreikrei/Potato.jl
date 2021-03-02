# =========================================================================
#    BASIC DATA GENERATION
# =========================================================================

const vertex_data = Ref{Any}(nothing)
V() = sort!(collect(keys(vertex_data[])))
V(i) = vertex_data[][i]

const distance_data = Ref{Any}(nothing)
dist() = distance_data[]
dist(i,j) = distance_data[][i,j]

const vehicle_data = Ref{Any}(nothing)
K() = sort!(collect(keys(vehicle_data[])))
K(k) = vehicle_data[][k]

const demand_data = Ref{Any}(nothing)
d() = demand_data[]
d(i,t) = demand_data[][i,t]

function extract!(path::String;t::Int64,f::String) #extract from excel
    xf = XLSX.readxlsx(path) #READ WORKSHEET
    data = Dict{Symbol,DataFrame}() #DATAFRAME DICT

    for sheets in XLSX.sheetnames(xf) #TURN SHEETS INTO DATAFRAME
        df = DataFrame(XLSX.gettable(xf[sheets])...) #TRANSFORM TABLE INTO DATAFRAME
        data[Symbol(sheets)] = df #DEFINE THE NAME FROM THE WORKSHEET
    end

    V = Dict{Int64,vtx}() #INITIATE VERTICES
    for v in eachrow(data[:vertices]) #ITERATE OVER DATA
        V[v.id] = vtx(
            v.name, v.type,
            v.x, v.y, v.MAX, v.MIN, v.START,
            v.h
        )
    end

    dist = JuMP.Containers.DenseAxisArray{Float64}(undef, keys(V), keys(V))
    dist .= 999999999
    for i in keys(V), j in keys(V)
        if i != j
            if f == "haversine"
                dist[i,j] = haversine([V[i].x,V[i].y],[V[j].x,V[j].y],6378.137)
            else #if f == "euclidean"
                dist[i,j] = euclidean([V[i].x,V[i].y],[V[j].x,V[j].y])
            end
        end
    end

    K = Dict{Int64,veh}() #INITIATE VEHICLES
    for k in eachrow(data[:vehicles]) #ITERATE OVER DATA
        K[k.id] = veh(
            k.name, k.type,
            parse.(Int64,split(k.cover)),
            Dict(parse.(Int64,split(k.cover)) .=>  parse.(Int64,split(k.BP))), k.Q,
            k.vx, k.vl, k.fp, k.fd
        )
    end

    d = JuMP.Containers.DenseAxisArray(
        Array{Float64}(data[:demands][:,string.(t)]), #dataset
        Array{Int64}(data[:demands].point)
    ) #one dimension

    vertex_data[] = V
    distance_data[] = dist
    vehicle_data[] = K
    demand_data[] = d

    return V,dist,K,d
end

function initStab()
    slackCoeff = sl_C()
    surpCoeff = su_C()
    slackLim = abs.(d())
    surpLim = abs.(d())

    return stabilizer(slackCoeff,surpCoeff,slackLim,surpLim)
end

function root()
    root = node(
        uuid1(), uuid1(),
        Vector{bound}(),Vector{col}(),
        initStab(),
        ["UNVISITED"]
    )

    return root
end
