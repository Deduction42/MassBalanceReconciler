using Accessors
using LinearAlgebra
include("_PlantInfo.jl")


#=============================================================================
Abstract Interface ("value" and "stdev" must be fields)
=============================================================================#
abstract type AbstractMeas{S, T} end

eltype(::Type{<:AbstractMeas{S,T}}) where {S,T} = T

stateindex(m::AbstractMeas) = collect(m.stream)
stateindex(v::AbstractVector{<:Species})  = reduce(vcat, v)
stateindex(v::AbstractVector{<:Reaction}) = stoich_extent.(v)
standarderr(x::AbstractVector, m::AbstractMeas) = innovation(x,m)/m.stdev

meastype(::Type{M}) where M <: AbstractMeas = Base.typename(M).wrapper
meastype(m::AbstractMeas) = meastype(typeof(m))

setvalue(m::AbstractMeas, v) = @set m.value = v
getvalue(m::AbstractMeas) = m.value

getvalue(d::Dict, v) = d[v]
getvalue(d::Dict, v::AbstractVector) = map(Base.Fix1(getindex, d), v)
getvalue(d::Dict, v::Species{S}) where S = Species{S}(getvalue(d, v[:]))

readvalue(m::AbstractMeas{S,T}, d::Dict) where {S,T} = setvalue(m, getvalue(d, m.tag))
updatethermo(m::AbstractMeas, d::Dict{Symbol,<:ThermoState}) = m

function negloglik(x::AbstractVector{T}, m::AbstractVector{<:AbstractMeas}) where T <: Real
    RT = promote_type(T,Float64)
    if isempty(m)
        return zero(RT)
    else
        return sum(Base.Fix1(negloglik, x), m)
    end
end


#=============================================================================
Univariate measuremnts
=============================================================================#
abstract type AbstractSingleMeas{S,T} <: AbstractMeas{S,T} end

function negloglik(x::AbstractVector{T}, m::AbstractSingleMeas) where T <: Real
    return 0.5*(innovation(x, m)/m.stdev)^2
end

function innovation(x::AbstractVector, m::AbstractSingleMeas)
    return m.value - prediction(x, m)
end

function innovation(x::AbstractVector{T}, vm::AbstractVector{AbstractSingleMeas}) where T <: Real
    return map(Base.Fix1(innovation, x), vm)
end

noisecov(m::AbstractSingleMeas) = m.stdev^2

#=============================================================================
Multivariate measuremnts
=============================================================================#
abstract type AbstractMultiMeas{S,T} <: AbstractMeas{S,T} end

function negloglik(x::AbstractVector{T}, m::AbstractMultiMeas) where T <: Real
    return 0.5*sum(x->x*x, innovation(x, m)./m.stdev)
end

function innovation(x::AbstractVector, m::AbstractMultiMeas)
    return m.value[:] .- prediction(x, m)
end

function innovation(x::AbstractVector{T}, vm::AbstractVector{AbstractMultiMeas{S}}) where {S, T<:Real}
    result = promote_type(T, Float64)[]
    for m in vm
        append!(result, innovation(x,m))
    end
    return result
end

noisecov(m::AbstractMultiMeas) = Diagonal(m.stdev[:].^2)

#=============================================================================
Volumetric flow rates
=============================================================================#
@kwdef struct VolumeFlowMeas{S, T, N} <: AbstractSingleMeas{S, T}
    id       :: Symbol
    streamid :: Symbol
    tag      :: String
    value    :: T
    molarvol :: Species{S, Float64, N}
    stream   :: Species{S, Int, N}
    stdev    :: Float64
end
VolumeFlowMeas{S, T}(x...) where {S,T} = VolumeFlowMeas{S, T, length(S)}(x...)
VolumeFlowMeas{S, T}(;kw...) where {S,T} = VolumeFlowMeas{S, T, length(S)}(;kw...)

function prediction(x::AbstractVector, m::VolumeFlowMeas)
    stream = x[m.stream[:]]
    return dot(m.molarvol[:], stream)
end

function build(::Type{<:VolumeFlowMeas}, measinfo::MeasInfo, streams::Dict{Symbol, <:StreamInfo{S}}, thermo::Dict{Symbol, <:ThermoState}) where S
    if length(measinfo.tags) != 1
        error("Measurement Type: VolumeFlowMeas only supports 1 tag, measurement id '$(measinfo.id)' contains $(length(measinfo.tags))")
    end

    thermostate = thermo[measinfo.stream]
    return VolumeFlowMeas{S, Float64}(
        id       = measinfo.id,
        streamid = measinfo.stream,
        tag      = measinfo.tags[1],
        value    = 0.0,
        stdev    = measinfo.stdev[1],
        stream   = streams[measinfo.stream].index,
        molarvol = molar_volumes(Species{S}, thermostate)
    )
end

function updatethermo(meas::VolumeFlowMeas{S}, thermo::Dict{Symbol, <:ThermoState}) where S
    thermostate = thermo[meas.streamid]
    return @set meas.molarvol = molar_volumes(Species{S}, thermostate)
end

#=============================================================================
Mass flow rates
=============================================================================#
@kwdef struct MassFlowMeas{S, T, N} <: AbstractSingleMeas{S, T}
    id          :: Symbol
    streamid    :: Symbol
    tag         :: String
    value       :: T
    molarmass   :: Species{S, Float64, N}
    stream      :: Species{S, Int, N}
    stdev       :: Float64
end
MassFlowMeas{S, T}(x...) where {S,T}  = MassFlowMeas{S, T, length(S)}(x...)
MassFlowMeas{S, T}(;kw...) where {S,T}  = MassFlowMeas{S, T, length(S)}(;kw...)

function prediction(x::AbstractVector, m::MassFlowMeas)
    stream = x[m.stream[:]]
    return dot(m.molarmass[:], stream)
end

function build(::Type{<:MassFlowMeas}, measinfo::MeasInfo, streams::Dict{Symbol, <:StreamInfo{S}}, thermo::Dict{Symbol,<:ThermoState}) where S
    if length(measinfo.tags) != 1
        error("Measurement Type: MassFlowMeas only supports 1 tag, measurement id '$(measinfo.id)' contains $(length(measinfo.tags))")
    end
    
    thermostate = thermo[measinfo.stream]
    return MassFlowMeas{S, Float64}(
        id        = measinfo.id,
        streamid  = measinfo.stream,
        tag       = measinfo.tags[1],
        value     = 0.0,
        stdev     = measinfo.stdev[1],
        stream    = streams[measinfo.stream].index,
        molarmass = molar_weights(Species{S}, thermostate)
    )
end

function updatethermo(meas::MassFlowMeas{S}, thermo::Dict{Symbol, <:ThermoState}) where S
    thermostate = thermo[meas.streamid]
    return @set meas.molarmass = molar_weights(Species{S}, thermostate)
end

#=============================================================================
Molar Analysis
=============================================================================#
@kwdef struct MoleAnalyzer{S, T, N} <: AbstractMultiMeas{S, T}
    id       :: Symbol
    streamid :: Symbol
    tag      :: Species{S, String, N}
    value    :: Species{S, T, N}
    stream   :: Species{S, Int, N}
    stdev    :: Species{S, Float64, N}
end
MoleAnalyzer{S, T}(x...) where {S,T} = MoleAnalyzer{S, T, length(S)}(x...)
MoleAnalyzer{S, T}(;kw...) where {S,T} = MoleAnalyzer{S, T, length(S)}(;kw...)

function prediction(x::AbstractVector, m::MoleAnalyzer)
    stream = x[m.stream[:]]
    return stream./sum(stream)
end

function build(::Type{<:MoleAnalyzer}, measinfo::MeasInfo, streams::Dict{Symbol, <:StreamInfo{S}}, thermo::Dict{Symbol,<:ThermoState}) where S
    measid = measinfo.id
    N = length(S)

    if length(measinfo.tags) != N
        error("Measurement Type: MassFlowMeas only supports $(N) tags, measurement id '$(measid)' contains $(length(measinfo.tags))")
    end

    return MoleAnalyzer{S, Float64}(
        id       = measid,
        streamid = measinfo.stream,
        tag      = Species{S,String}(measinfo.tags),
        value    = zero(Species{S,Float64,N}),
        stdev    = Species{S,Float64,N}(measinfo.stdev),
        stream   = streams[measinfo.stream].index
    )
end


#=============================================================================
Mole Balancer
=============================================================================#
@kwdef struct MoleBalance{S, T, N} <: AbstractMultiMeas{S, T}
    id        :: Symbol
    value     :: Species{S, T, N}
    interval  :: Float64
    inlets    :: Vector{Species{S, Int, N}}
    outlets   :: Vector{Species{S, Int, N}}
    reactions :: Vector{Reaction{S, Int, N}}
    stdev     :: Species{S, Float64, N}
end
MoleBalance{S, T}(x...) where {S,T} = MoleBalance{S, T, length(S)}(x...)
MoleBalance{S, T}(;kw...) where {S,T} = MoleBalance{S, T, length(S)}(;kw...)

function stateindex(m::MoleBalance)
    return [
        stateindex(m.inlets);
        stateindex(m.outlets);
        stateindex(m.reactions)
    ]
end

function prediction(x::AbstractVector{T}, m::MoleBalance{S, <:Float64, N}) where {S,T,N}
    RT = promote_type(T,Float64)
    default = zero(SVector{N,RT})

    balance = (
          (isempty(m.inlets)    ? default : sum(Base.Fix1(speciesvec, x), m.inlets))
        - (isempty(m.outlets)   ? default : sum(Base.Fix1(speciesvec, x), m.outlets))
        + (isempty(m.reactions) ? default : sum(Base.Fix1(speciesvec, x), m.reactions))
    )

    return balance.*m.interval
end

function build(::Type{<:MoleBalance}, nodeinfo::NodeInfo, streams::Dict{Symbol, <:StreamInfo{S}}) where S
    nodeid = nodeinfo.id
    N = length(S)

    return MoleBalance{S, Float64}(
        id        = nodeid,
        value     = zero(Species{S, Float64, N}),
        interval  = 0.0,
        stdev     = Species{S}(nodeinfo.stdev),
        inlets    = [streams[id].index for id in nodeinfo.inlets],
        outlets   = [streams[id].index for id in nodeinfo.outlets],
        reactions = nodeinfo.reactions
    )
end

function setinterval(m::MoleBalance{S,T}, Δt::Real) where {S,T} 
    return @set m.interval = T(Δt)
end

function setintervals!(vm::AbstractVector{<:MoleBalance}, Δt::Real)
    vm .= setinterval.(vm, Δt)
    return vm
end

#=============================================================================
Collection of all measurements
=============================================================================#
@kwdef struct MeasCollection{S,T,N}
    VolumeFlowMeas  :: Vector{VolumeFlowMeas{S,T,N}} = VolumeFlowMeas{S,T,N}[]
    MassFlowMeas    :: Vector{MassFlowMeas{S,T,N}}   = MassFlowMeas{S,T,N}[]
    MoleAnalyzer    :: Vector{MoleAnalyzer{S,T,N}}   = MoleAnalyzer{S,T,N}[]
    MoleBalance     :: Vector{MoleBalance{S,T,N}}    = MoleBalance{S,T,N}[]
end
MeasCollection{S,T}(;kwargs...) where {S,T}= MeasCollection{S,T,length(S)}(kwargs...)

Base.getindex(m::MeasCollection, k::Symbol)  = getproperty(m, k)
Base.getindex(m::MeasCollection, k::Integer) = getproperty(m, fieldnames(MeasCollection)[k])
Base.getindex(m::MeasCollection, k::Colon)   = map(Base.Fix1(getproperty, m), fieldnames(MeasCollection))
Base.firstindex(m::MeasCollection) = 1
Base.lastindex(m::MeasCollection) = length(fieldnames(MeasCollection))
Base.getindex(m::MeasCollection, ::Type{T}) where T = getproperty(m, Symbol(T))
Base.getindex(m::MeasCollection, k::AbstractVector) = map(Base.Fix1(getindex, m), k)
Base.getindex(m::MeasCollection, k::Tuple) = map(Base.Fix1(getindex, m), k)

function readvalues!(vmeas::AbstractVector{M}, d::Dict) where {M <: AbstractMeas}
    reader = Base.Fix2(readvalue, d)
    if hasfield(M, :tag)
        vmeas .= reader.(vmeas)
    end
    return vmeas
end

function readvalues!(c::MeasCollection, d::Dict)
    for fn in fieldnames(MeasCollection)
        readvalues!(c[fn], d)
    end
    return c
end

function updatethermo!(vmeas::AbstractVector{M}, d::Dict{Symbol, <:ThermoState}) where {M <: AbstractMeas}
    updater = Base.Fix2(updatethermo, d)
    vmeas .= updater.(vmeas)
    return vmeas
end

function updatethermo!(c::MeasCollection, d::Dict{Symbol, <:ThermoState})
    for fn in fieldnames(MeasCollection)
        updatethermo!(c[fn], d)
    end
    return c
end

function negloglik(x::AbstractVector, c::MeasCollection)
    return sum(fn-> negloglik(x, c[fn]), fieldnames(MeasCollection))
end
