include("_MeasCollection.jl")

using LinearAlgebra
using Dates
using TimeRecords
using OSQP
const TIMESTAMP_KEY = "_UNIX_TIMESTAMP"


@kwdef struct PlantClock
    timestamp :: Base.RefValue{Float64}
    interval  :: Base.RefValue{Float64}
    stepsize  :: Float64
end

@kwdef struct PlantState{L, N}
    clock        :: PlantClock
    thermo       :: ThermoModel{L,N}
    tags         :: Vector{TagInfo}
    statevec     :: Vector{Float64}
    stateinv     :: SparseMatrixCSC{Float64, Int64}
    predictor    :: @NamedTuple{A::Matrix{Float64}, iQ::Matrix{Float64}}
    measurements :: MeasCollection{L, Float64, N}
    streams      :: Vector{StreamRef{L,N}}
    nodes        :: Vector{NodeRef{L,N}}
    solver       :: OSQP.Model
end

@kwdef struct PlantSeries{L,N}
    plant  :: PlantState{L,N}
    states :: TimeSeries{Vector{Float64}}
    stdevs :: TimeSeries{Vector{Float64}}
end

PlantSeries(plant::PlantState{L,N}, states, stdevs) where {L,N} = PlantSeries{L,N}(plant, states, stdevs)

function PlantState(plantinfo::PlantInfo)
    #Retrieve the species vector (the main plant parameter)
    L = Tuple(plantinfo.thermo.labels)
    N = length(L)
    return PlantState{L,N}(plantinfo)
end

PlantState{L}(plantinfo::PlantInfo) where {L} = PlantState{L,length(L)}(plantinfo)

function PlantState{L,N}(plantinfo::PlantInfo) where {L,N}
    #Fill the composition copy slots (for shortcut constraints)
    fill_copycomps!(plantinfo)

    #Build the plant clock information
    stepsize = plantinfo.interval
    interval = Ref(plantinfo.interval)
    timestamp = Ref(datetime2unix(floor(now(UTC), Day(1))))

    plantclock = PlantClock(
        timestamp = timestamp,
        interval = interval,
        stepsize = stepsize,
    )

    #Build the initial state index
    indref   = Ref(0)

    #Build the main thermodynamic model
    thermo  = ThermoModel{L}(plantinfo.thermo)

    #Build the streams and index them
    streams = [StreamRef{L}(stream) for stream in plantinfo.streams]
    stateindex!(streams, indref)
    streamdict = Dict(stream.id=>stream for stream in streams)

    #Build the nodes and index them (mostly if they have reactions)
    nodes = [NodeRef{L}(node, streamdict) for node in plantinfo.nodes]
    stateindex!(nodes, indref)

    #The length of the state is the final index value
    Nx = indref[]

    #Initialize the state vector and make sure there are no zeros
    statevec = build_state(Nx, thermo, streams, plantinfo.streams)
    statevec = max.(statevec, 1e-9)

    #Build the state covariance assuming the nominal values are the standard deviation
    transmat   = Matrix(Diagonal(ones(Nx)))
    stateinv   = spzeros(Nx,Nx)
    predictor  = (A=transmat, iQ=copy(stateinv))

    #Build the measurements based off the thermodynamic information
    meascollection = MeasCollection{L,Float64}()

    for measinfo in plantinfo.measurements
        Mtype = measinfo.type
        meas = Mtype(measinfo, streamdict, thermo)
        push!(meascollection[Mtype], meas)
    end

    for nodeinfo in plantinfo.nodes
        meas = MoleBalance(nodeinfo, streamdict, interval)
        push!(meascollection[MoleBalance], meas)
    end

    assign_flowrefs!(meascollection)

    #Populate the final object with constructed values and pass through the stream and node information
    plantstate = PlantState{L,N}(
        clock = plantclock,
        thermo = thermo,
        tags = gettags(meascollection),
        statevec = statevec,
        stateinv = stateinv,
        predictor = predictor,
        measurements = meascollection,
        streams = streams,
        nodes = nodes,
        solver = OSQP.Model()
    )

    #Set up the solver problem
    solver_setup!(plantstate)

    return plantstate
end

tagnames(plant::PlantState) = getname.(plant.tags)

function getrecords(plantseries::PlantSeries{L}, stream::Symbol; components=L) where L
    streamref = getstreamref(plantseries, stream)
    return getrecords(plantseries, streamref, components=components)
end

function getrecords(plantseries::PlantSeries{L}, streamref::StreamRef; components=L) where L
    compinds  = map(c->streamref.index[c], components)
    return (
        time  = timestamps(plantseries.states),
        state = map(ind->[x.v[ind] for x in plantseries.states], compinds),
        stdev = map(ind->[σ.v[ind] for σ in plantseries.stdevs], compinds)
    )
end

function getrecords(plantseries::PlantSeries{L}, streamref::Nothing; components=L) where L
    return (
        time  = timestamps(plantseries.states),
        state = map(ind->zeros(length(plantseries.states)), components),
        stdev = map(ind->zeros(length(plantseries.states)), components)
    )
end

function getrecords(plantstate::PlantState{L}, stream::Symbol; components=L) where L
    streamref = getstreamref(plantstate, stream)
    return getrecords(plantstate, streamref, components=components)
end

function getrecords(plantstate::PlantState{L}, streamref::StreamRef; components=L) where L
    compinds  = map(c->streamref.index[c], components)
    return (
        time  = plantstate.clock.timestamp,
        state = plantstate.statevec[compinds],
        stdev = inv.(sqrt.(diag(plantstate.stateinv)[compinds]))
    )
end

function getrecords(plantstate::PlantState{L}, streamref::Nothing; components=L) where L
    return (time=plantstate.clock.timestamp, state=map(x->0.0, components), stdev=map(x->0.0, components))
end

getstreamref(plant::PlantSeries, id::Symbol) = getstreamref(plant.plant, id)

function getstreamref(plant::PlantState, id::Symbol)
    for s in plant.streams
        if s.id == id 
            return s
        end
    end
    throw(ArgumentError("No streams found in plant with id '$(id)'"))
end


function predict!(plant::PlantState)
    info_halflife = 300.0
    plant.stateinv .= plant.stateinv.*(2^(-plant.clock.interval[]/info_halflife))

    #Predict the state
    #A = plant.predictor.A
    #plant.statevec .= A*plant.statevec

    #Predict the state covariance using the Woodbury Formula
    #=
    #https://tlienart.github.io/posts/2018/12/13-matrix-inversion-lemmas/index.html
    #inv(F*G*H + E) ≈ iE - iE*F*inv(iG + H*iE*F)*H*iE
    iP = plant.stateinv
    iQ = Diagonal(plant.predictor.iQ)
    F = A 
    H = A'
    plant.stateinv .= iQ .- iQ*F*inv(hermitianpart(iP + H*iQ*F))*H*iQ
    hermitianpart!(plant.stateinv)
    return plant 
    =#
end

function loglik(x::AbstractVector, plant::PlantState)
    Δx = x .- plant.statevec
    state_loglik = -0.5*Δx'*plant.stateinv*Δx
    meas_loglik  = loglik(x, plant.measurements)
    return state_loglik + meas_loglik
end

function readvalues!(plant::PlantState, data::AbstractDict{<:AbstractString,<:Real})
    setclock!(plant, data)
    readvalues!(plant.measurements, data)
    calcvalues!(plant.measurements, plant.thermo)
    return plant
end

function updatethermo!(plant::PlantState)
    updatethermo!(plant.measurements, plant.statevec, plant.thermo)
    return plant 
end


function setclock!(plant::PlantState, data::AbstractDict{<:AbstractString}; nominal_interval=false)
    timestamp = data[TIMESTAMP_KEY]

    if nominal_interval
        plant.clock.interval[] = plant.clock.stepsize

    elseif (timestamp <= plant.clock.timestamp[])
        @warn "Dataset is later than the current plant state, assuming nominal interval"
        plant.clock.interval[] = plant.clock.stepsize

    else
        plant.clock.interval[] = (timestamp - plant.clock.timestamp[])
    end

    plant.clock.timestamp[] = timestamp
    return plant 
end


#==============================================================================================================================
Fill state vector with stream information defaults
==============================================================================================================================#
@kwdef struct WeightedValue{T}
    weight :: Float64
    value  :: T 
end

weight(x::WeightedValue) = x.weight
value(x::WeightedValue)  = x.value

function build_state(N::Integer, model::ThermoModel{L}, streamrefs::Vector{<:StreamRef{L}}, streaminfo::Vector{StreamInfo}) where L
    weighted_state = fill(WeightedValue(0.0, NaN), N)
    molweights = molar_weights(model)

    if length(streamrefs) != length(streaminfo)
        error("streamrefs and streaminfo must have same lengths")
    end

    #Fill all streams that have no parent
    for (streamref, streaminfo) in zip(streamrefs, streaminfo)
        if !hasparents(streamref)
            molefracs = Species{L}(fractions([streaminfo.molefracs[l] for l in L]))
            molweight = molefracs' * molweights.data
            moleflows = molefracs*(streaminfo.massflow/molweight)

            for (ind, flow) in zip(streamref.comp, moleflows)
                weighted_state[ind] = update(weighted_state[ind], flow)
            end
        end
    end

    #Fill all streams that have a parent
    for (streamref, streaminfo) in zip(streamrefs, streaminfo)
        if hasparents(streamref)
            moleflows = value.(speciesvec(weighted_state, streamref.comp))
            molefracs = fractions(moleflows)
            molweight = molefracs' * molweights.data
            totalflow = streaminfo.massflow/molweight

            ind = streamref.flow
            weighted_state[ind] = update(weighted_state[ind], totalflow)
        end
    end

    #Check for missing compositions and flows 
    state = value.(weighted_state)
    missing_comps = Symbol[]
    missing_flows = Symbol[]

    for streamref in streamrefs
        if hasparents(streamref)
            if isnan(state[streamref.flow])
                push!(missing_flows, streamref.id)
            end
        else
            if any(isnan, state[streamref.comp.data])
                push!(missing_comps, streamref.id)
            end
        end
    end

    if (!isempty(missing_comps)) & (!isempty(missing_flows))
        error("Incomplete system design info:\nThe following streams are missing compositions: $(missing_comps) \nThe following streams are missing mass flows: $(missing_flows)")
    elseif !isempty(missing_comps)
        error("Incomplete system design info:\nThe following streams are missing compositions: $(missing_comps)")
    elseif !isempty(missing_flows)
        error("Incomplete system design info:\nThe following streams are missing mass flows: $(missing_flows)")
    end

    #Any remaining states are set to zero
    state[isnan.(state)] .= 0

    return state
end


function update(oldval::WeightedValue{T}, newval) where T
    if iszero(oldval.weight)
        return WeightedValue{T}(weight=1, value=newval)
    end

    newweight = 1 + oldval.weight
    return WeightedValue{T}(
        weight = newweight, 
        value = (oldval.weight*oldval.value + newval)/newweight
    )
end


#==============================================================================================================================
Set state/transition covariance by simulating data and estimating a maximum ignorance problem
==============================================================================================================================#
#=
function set_covariances!(plant::PlantState)
    X = plant.statevec
    meascol  = plant.measurements
    tagpairs = Pair{String,Float64}[]

    #Fill out the measurements
    for fn in fieldnames(MeasCollection)
        if fn != :MoleBalance
            for meas in meascol[fn]
                append!(tagpairs, simulate_measurements(X, meas))
            end
        end
    end

    #Fill out the expected timestamp
    unix_t = (plant.clock.timestamp[] + plant.clock.interval[])
    push!(tagpairs, TIMESTAMP_KEY => unix_t)

    #Fill in the measurement values
    meas_dict = Dict(tagpairs)
    readvalues!(plant, meas_dict)
    updatethermo!(plant)

    #Calculate the state covariance at this (close to optimal) location
    P⁻¹ = reconcile_statecov!(plant)

    #Use this posterior covariance to scale the initial state transition matrix
    iQ = Matrix(0.5*Diagonal(P⁻¹))
    plant.predictor.iQ .= iQ

    return plant
end

function simulate_measurements(X::AbstractVector, meas::AbstractSingleMeas)
    return SVector(meas.tag => prediction(X, meas))
end

function simulate_measurements(X::AbstractVector, meas::AbstractMultiMeas)
    return speciesvec(meas.tag) .=> prediction(X, meas)
end

function simulate_measurements(X::AbstractVector, meas::VolumeFlowMeas)
    return SVector(
        meas.tag[:V] => prediction(X, meas),
        string(meas.tag[:T]) => meas.value[:T],
        string(meas.tag[:P]) => meas.value[:P]
    )
end
=#