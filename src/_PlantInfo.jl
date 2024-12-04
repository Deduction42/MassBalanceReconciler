include("_ThermoModel.jl")
using Accessors


function stateindex!(indref::Base.RefValue, s::Species{L}) where {L}
    N = length(L)
    start = indref[] + 1
    indref[] = indref[] + N
    return Species{L}(SVector{N}(start:indref[]))
end

function stateindex!(indref::Base.RefValue, r::Integer)
    indref[] = indref[] + 1
    return indref[]
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
    massflow :: Float64
    index :: Species{L, Int, N}
    phase :: Symbol = :unknown
    refid :: Symbol = :nothing
    scale :: Int = 0
end

function StreamInfo{L}(;id, massflow, refid=:nothing, phase=:unknown) where L
    N = length(L)
    return StreamInfo{L,N}(
        id = id,
        index = Species{L, Int, N}(zero(SVector{N,Int})),
        massflow = massflow,
        phase = phase
    )
end

function stateindex!(indref::Base.RefValue, streaminfo::StreamInfo{L}) where {L}
    if streaminfo.refid == :nothing
        return @set streaminfo.index = stateindex!(indref, streaminfo.index)
    else
        return @set streaminfo.scale = stateindex!(indref, streaminfo.scale)
    end
end

function Base.getindex(v::AbstractVector, ind::StreamInfo)
    streamvec = v[ind.index[:]]
    if ind.refid == :nothing
        return streamvec
    else
        return streamvec.*v[ind.scale]
    end
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
    stream :: Symbol = :nothing
    node   :: Symbol = :nothing
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

    _fillrefs!(plantinfo.streams)

    return indref[]
end

function _fillrefs!(streams::Vector{<:StreamInfo})
    streamdict = Dict( stream.id=>stream for stream in streams)

    for (ii, stream) in enumerate(streams)
        if stream.refid != :nothing
            streams[ii] = @set stream.index = streamdict[stream.refid].index 
        end
    end
    return streams
end



#=============================================================================
Thermodynamic information (separated from plant to enable abstraction of composition)
=============================================================================#
@kwdef struct ThermoInfo{L,N}
    tags   :: Dict{Symbol, ThermoState{L,String,N}}
    values :: Dict{Symbol, ThermoState{L,Float64,N}}
end
ThermoInfo{L}(x...) where L = ThermoInfo{L,length(L)}(x...)
ThermoInfo{L}(;kw...) where L = ThermoInfo{L,length(L)}(;kw...)

function readvalues!(obj::ThermoInfo, d::Dict)
    for (k, tags) in obj.tags
        obj.values[k] = readvalues(d, tags)
    end
    return obj
end

