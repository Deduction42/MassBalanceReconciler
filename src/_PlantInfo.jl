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
    phase :: Symbol = :unknown
end

function StreamInfo{L}(;id, phase=unknown) where L
    return StreamInfo{L,N}(
        id = id,
        index = Species{L, Int, N}(zero(SVector{N,Int})),
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
Construction info for entire system
=============================================================================#
@kwdef struct PlantInfo{L,N}
    streams :: Dict{Symbol, StreamInfo{L,N}} = Dict{Symbol, StreamInfo{L,N}}()
    nodes   :: Dict{Symbol, NodeInfo{L,N}}   = Dict{Symbol, NodeInfo{L,N}}()
    measurements :: Dict{Symbol, MeasInfo}   = Dict{Symbol, NodeInfo{L,N}}()
end

PlantInfo{L}(;kwargs...) where {L} = PlantInfo{L, length(L)}(;kwargs...)


function stateindex!(plantinfo::PlantInfo)
    indref = Ref(0)

    for k in sort!(collect(keys(plantinfo.streams)))
        plantinfo.streams[k] = stateindex!(indref, plantinfo.streams[k])
    end

    for k in sort!(collect(keys(plantinfo.nodes)))
        plantinfo.nodes[k] = stateindex!(indref, plantinfo.nodes[k])
    end

    return indref[]
end