include("_AbstractMeas.jl")










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
    species :: Vector{Symbol}
    thermo  :: Dict{Symbol,Dict{String,Float64}}
    streams :: Vector{StreamInfo} = StreamInfo[]
    nodes   :: Vector{NodeInfo}   = NodeInfo[]
    measurements  :: Vector{MeasInfo}  = MeasInfo[]
    relationships :: Vector{StreamRelationship} = StreamRelationship[]
end

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

