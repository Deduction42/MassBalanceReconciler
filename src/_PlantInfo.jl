include("_ThermoModel.jl")
using Accessors

function stateindex!(indref::Base.RefValue, s::Species{L}) where {L}
    N = length(L)
    start = indref[] + 1
    indref[] = indref[] + N
    return Species{L}(SVector{N}(start:indref[]))
end

function stateindex!(indref::Base.RefValue, r::Reaction{L}) where {L}
    indref[] = indref[] + 1
    return @set r.extent = indref[]
end

#=============================================================================
Construction info for streams
=============================================================================#
@kwdef struct StreamInfo{L, N}
    id :: Symbol
    index :: Species{L, Int, N}
    massflow :: Float64
    phase :: Symbol = :unknown
end

function StreamInfo{L}(;id, massflow, phase=unknown) where L
    return StreamInfo{L,N}(
        id = id,
        index = Species{L, Int, N}(zero(SVector{N,Int})),
        massflow = massflow,
        phase = phase
    )
end

function stateindex!(indref::Base.RefValue, streaminfo::StreamInfo{L}) where {L}
    return @set streaminfo.index = stateindex!(indref, streaminfo.index)
end

#=============================================================================
Construction info for nodes
=============================================================================#
@kwdef struct NodeInfo{L, N}
    id        :: Symbol
    stdev     :: Species{L,Float64,N}
    inlets    :: Vector{Symbol}
    outlets   :: Vector{Symbol}
    reactions :: Vector{Reaction{L, Int, N}}
end

function NodeInfo{L}(;id, inlets, outlets, stdev) where {L}
    N = length(L)
    return NodeInfo{L, N}(
        id = id,
        stdev  = stdev,
        inlets = inlets,
        outlets = outlets,
        reactions = Reaction{L, Int, N}[],
    )
end

function stateindex!(indref::Base.RefValue, nodeinfo::NodeInfo{L}) where {L}
    for (ii, reaction) in enumerate(nodeinfo.reactions)
        nodeinfo.reactions[ii] = stateindex!(indref, reaction)
    end
    return nodeinfo
end


function add_reaction!(nodeinfo::NodeInfo{L,N}, stoich::Species) where {L,N}
    push!(nodeinfo.reactions, Reaction{L,Int,N}(0, stoich))
end

#=============================================================================
Construction info for measurements
=============================================================================#
@kwdef struct MeasInfo
    id     :: Symbol
    type   :: UnionAll
    tags   :: Vector{String}
    stdev  :: Vector{Float64}
    stream :: Symbol = Symbol("")
    node   :: Symbol = Symbol("")
end

#=============================================================================
Construction info for simple stream relationshps
=============================================================================#
@kwdef struct StreamRelationship
    id     :: Symbol
    parent :: Symbol
    factor :: Float64
    timeconst :: Float64
end

#=============================================================================
Construction info for entire system
=============================================================================#
@kwdef struct PlantInfo{L,N}
    streams :: Vector{StreamInfo{L,N}} = StreamInfo{L,N}[]
    nodes   :: Vector{NodeInfo{L,N}}   = NodeInfo{L,N}[]
    measurements  :: Vector{MeasInfo}  = MeasInfo[]
    relationships :: Vector{StreamRelationship} = StreamRelationship[]
end

PlantInfo{L}(;kwargs...) where {L} = PlantInfo{L, length(L)}(;kwargs...)


function stateindex!(plantinfo::PlantInfo)
    indref = Ref(0)

    for k in eachindex(plantinfo.streams)
        plantinfo.streams[k] = stateindex!(indref, plantinfo.streams[k])
    end

    for k in eachindex(plantinfo.nodes)
        plantinfo.nodes[k] = stateindex!(indref, plantinfo.nodes[k])
    end

    return indref[]
end

#=============================================================================
Thermodynamic information (separated from plant to enable abstraction of composition)
=============================================================================#
@kwdef struct ThermoInfo{L,N}
    tags   :: Dict{Symbol, ThermoState{L,String,N}}
    values :: Dict{Symbol, ThermoState{L,Float64,N}}
end

function readvalues!(d::Dict, obj::ThermoInfo)
    for (k, tags) in obj.tags
        obj.values[k] = readvalues(d, tags)
    end
    return obj
end

