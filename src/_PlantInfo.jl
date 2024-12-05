include("_AbstractStreamRef.jl")


#=============================================================================
Construction info for streams
=============================================================================#
@kwdef struct StreamInfo
    id        :: Symbol
    massflow  :: Float64
    molefracs :: Union{Symbol, Dict{Symbol, Float64}}
end

#=============================================================================
Construction info for nodes
=============================================================================#
@kwdef struct NodeInfo
    id        :: Symbol
    stdev     :: Dict{Symbol, Float64}
    inlets    :: Vector{Symbol}
    outlets   :: Vector{Symbol}
    reactions :: Vector{Dict{Symbol, Float64}} = Dict{Symbol, Float64}[]
end

function add_reaction!(nodeinfo::NodeInfo{L,N}, stoich::Species) where {L,N}
    push!(nodeinfo.reactions, ReactionRef{L,N}(0, stoich))
end

#=============================================================================
Construction info for measurements
=============================================================================#
@kwdef struct MeasInfo
    id     :: Symbol
    type   :: UnionAll
    tags   :: Dict{Symbol, Union{String,Float64}}
    stdev  :: Dict{Symbol, Float64}
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
@kwdef struct PlantInfo
    streams :: Vector{StreamInfo} = StreamInfo[]
    nodes   :: Vector{NodeInfo}   = NodeInfo[]
    measurements  :: Vector{MeasInfo}  = MeasInfo[]
    relationships :: Vector{StreamRelationship} = StreamRelationship[]
end

PlantInfo{L}(;kwargs...) where {L} = PlantInfo{L, length(L)}(;kwargs...)

#=
function stateindex!(plantinfo::PlantInfo)
    indref = Ref(0)

    for k in eachindex(plantinfo.streams)
        plantinfo.streams[k] = stateindex!(indref, plantinfo.streams[k])
    end

    for k in eachindex(plantinfo.nodes)
        plantinfo.nodes[k] = stateindex!(indref, plantinfo.nodes[k])
    end

    fillrefs!(plantinfo.streams)

    return indref[]
end
=#





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

