include("_AbstractMeas.jl")
using LinearAlgebra
using Dates
using TimeRecords
const TIMESTAMP_KEY = "_UNIX_TIMESTAMP"

#=============================================================================
Construction info for entire system
=============================================================================#
@kwdef struct PlantInfo <: AbstractInfo
    interval :: Float64
    thermo  :: ThermoInfo
    streams :: Vector{StreamInfo} = StreamInfo[]
    nodes   :: Vector{NodeInfo}   = NodeInfo[]
    measurements  :: Vector{MeasInfo}  = MeasInfo[]
    relationships :: Vector{StreamRelationship} = StreamRelationship[]
end

function PlantInfo(d::AbstractDict{<:Symbol})
    return PlantInfo(
        interval = d[:interval],
        thermo  = ThermoInfo(d[:thermo]),
        streams = StreamInfo.(d[:streams]),
        nodes   = NodeInfo.(d[:nodes]),
        measurements = MeasInfo.(d[:measurements]),
        relationships = StreamRelationship.(d[:relationships])
    )
end

@kwdef struct PlantClock
    timestamp :: Base.RefValue{Float64}
    interval  :: Base.RefValue{Float64}
    stepsize  :: Float64
end

@kwdef struct PlantState{L, N}
    clock        :: PlantClock
    thermo       :: ThermoModel{L,N}
    statevec     :: Vector{Float64}
    statecov     :: Matrix{Float64}
    dpredictor   :: @NamedTuple{A::Matrix{Float64}, Q::Matrix{Float64}}
    measurements :: MeasCollection{L, Float64, N}
    streams      :: Vector{StreamRef{L,N}}
    nodes        :: Vector{NodeRef{L,N}}
end

@kwdef struct PlantSeries{L,N}
    plant  :: PlantState{L,N}
    states :: TimeSeries{Vector{Float64}}
    stdevs :: TimeSeries{Vector{Float64}}
end

PlantSeries(plant::PlantState{L,N}, states, stdevs) where {L,N} = PlantSeries{L,N}(plant, states, stdevs)
PlantSeries(;kwargs...) = PlantSeries(;kwargs...)

function PlantState(plantinfo::PlantInfo)

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

    #Retrieve the species vector (the main plant parameter)
    L = Tuple(plantinfo.thermo.labels)

    #Build the main thermodynamic model
    thermo  = ThermoModel{L}(plantinfo.thermo)

    #Build the streams and index them
    streams = [StreamRef{L}(stream) for stream in plantinfo.streams]
    stateindex!(indref, streams)
    streamdict = Dict(stream.id=>stream for stream in streams)

    #Build the nodes and index them (mostly if they have reactions)
    nodes = [NodeRef{L}(node, streamdict) for node in plantinfo.nodes]
    stateindex!(indref, nodes)

    #The length of the state is the final index value
    Nx = indref[]

    #Initialize the state vector
    statevec = zeros(Float64, Nx)

    #Fill the state vector according to stream flow defaults
    fillstate!(statevec, thermo, streams, plantinfo.streams)

    #Fill the state vector according to reaction stoichiometry
    fillstate!(statevec, nodes)

    #Build the state covariance assuming the nominal values are the standard deviation    #Obtain the state transition object
    transmat = state_transition(Nx, plantinfo.relationships, streamdict)
    statecov = Matrix(Diagonal(statevec.^2))
    dpredictor = (A=transmat, Q=statecov)


    #Build the measurements based off the thermodynamic information
    meascollection = MeasCollection{L,Float64}()

    for measinfo in plantinfo.measurements
        meastype = measinfo.type
        meas = meastype(measinfo, streamdict, thermo)
        push!(meascollection[meastype], meas)
    end

    for nodeinfo in plantinfo.nodes
        meas = MoleBalance(nodeinfo, streamdict, interval)
        push!(meascollection[MoleBalance], meas)
    end

    #Populate the final object with constructed values and pass through the stream and node information
    return PlantState{L,length(L)}(
        clock = plantclock,
        thermo = thermo,
        statevec = statevec,
        statecov = statecov,
        dpredictor = dpredictor,
        measurements = meascollection,
        streams = streams,
        nodes = nodes
    )
end


function predict!(plant::PlantState)
    interval = plant.clock.interval[]
    Ad = exp(interval.*plant.dpredictor.A)
    plant.statevec .= Ad*plant.statevec
    plant.statecov .= Ad'*plant.statecov*Ad + interval.*plant.dpredictor.Q
    return plant 
end

function negloglik(x::AbstractVector, plant::PlantState)
    Δx = x .- plant.statevec
    state_negloglik = Δx'*plant.statecov*Δx
    meas_negloglik  = negloglik(x, plant.measurements)
    return state_negloglik + meas_negloglik
end

function readvalues!(plant::PlantState, data::AbstractDict{<:AbstractString,<:Real})
    setclock!(plant, data)
    readvalues!(plant.measurements, data)
    return plant
end

function updatethermo!(plant::PlantState)
    updatethermo!(plant.measurements, plant.statevec, plant.thermo)
    return plant 
end

#==============================================================================================================================
Fill state vector with stream information defaults
==============================================================================================================================#

function fillstate!(X::AbstractVector, model::ThermoModel, streamrefs::Vector{<:StreamRef}, streaminfo::Vector{StreamInfo})
    molweights = molar_weights(model)

    if length(streamrefs) != length(streaminfo)
        error("streamrefs and streaminfo must have same lengths")
    end

    #Fill all streams that have no parent
    for (streamref, streaminfo) in zip(streamrefs, streaminfo)
        if !hasparent(streamref)
            fillstate!(X, molweights, streamref, streaminfo)
        end
    end

    #Fill all streams that have a parent
    for (streamref, streaminfo) in zip(streamrefs, streaminfo)
        if hasparent(streamref)
            fillstate!(X, molweights, streamref, streaminfo)
        end
    end

    return X
end


function fillstate!(X::AbstractVector, molweights::Species{L}, streamref::StreamRef{L}, streaminfo::StreamInfo) where L
    if streamref.id != streaminfo.id
        error("Identifiers of streaminfo ($(streaminfo.id)) and streamref ($(streamref.id)) must match")
    end

    if hasparent(streaminfo)
        refmols = speciesvec(X, streamref.index)
        refmass = dot(molweights, refmols)
        X[streamref.scale] = streaminfo.massflow/refmass
    else
        molefracs = [streaminfo.molefracs[l] for l in L]
        molemass  = dot(molweights, molefracs./sum(molefracs))
        X[streamref.index[:]] = molefracs.*(streaminfo.massflow/molemass)
    end

    return X
end

function fillstate!(X::AbstractVector, noderefs::AbstractVector{<:NodeRef{L}}) where L
    for noderef in noderefs
        fillstate!(X, noderef)
    end
    return X 
end

function fillstate!(X::AbstractVector, noderef::NodeRef{L}) where L
    inputs = Species{L}(sum(Base.Fix1(speciesvec, X), noderef.inlets))

    for reaction in noderef.reactions
        X[reaction.extent] = stoich_extent(reaction, inputs)
    end
    return X 
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











