using Accessors
using LinearAlgebra
include("_AbstractStreamRef.jl")



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

function MeasInfo(d::AbstractDict{Symbol}) 
    type = eval(Meta.parse(d[:type]))
    return MeasInfo(
        id     = Symbol(d[:id]),
        type   = type,
        tags   = symbolize(Union{String,Float64}, d[:tags]),
        stdev  = symbolize(Float64, d[:stdev]),
        stream = Symbol(get(d, :stream, :nothing)),
        node   = Symbol(get(d, :node, :nothing))
    )
end

MeasInfo(d::AbstractDict{<:AbstractString}) = MeasInfo(symbolize(d))

#=
function build(info::MeasInfo, streams::Dict{Symbol, StreamRef})
    return build(info.type, info, streams)
end
=#

#=============================================================================
Abstract Interface ("value" and "stdev" must be fields)
=============================================================================#
abstract type AbstractMeas{S, T} end

eltype(::Type{<:AbstractMeas{S,T}}) where {S,T} = T
meastype(::Type{M}) where M <: AbstractMeas = Base.typename(M).wrapper
meastype(m::AbstractMeas) = meastype(typeof(m))

#state index collection =======================================================
function addinds!(inds::BitArray, v::AbstractVector{<:Integer})
    inds[v] .= true
    return inds 
end

function addinds!(inds::BitArray, ind::Integer)
    inds[ind] = true
    return inds 
end

function addinds!(inds::BitArray, s::StreamRef) 
    addinds!(inds, s.index[:])
    if s.refid != :nothing
        addinds!(inds, s.scale)
    end
    return inds
end

function addinds!(inds::BitArray, v::AbstractVector)
    for x in v
        addinds!(inds, x)
    end
    return inds 
end

addinds!(inds::BitArray, r::ReactionRef{L}) where L = addinds!(inds, r.extent)
addinds!(inds::BitArray, m::AbstractMeas) = addinds!(inds, m.stream)


#tag-reading interface ===========================================================
setvalue(m::AbstractMeas, v) = @set m.value = v
getvalue(m::AbstractMeas) = m.value

getvalue(d::Dict, v) = d[v]
getvalue(d::Dict, v::AbstractVector) = map(Base.Fix1(getindex, d), v)
getvalue(d::Dict, v::Species{S}) where S = Species{S}(getvalue(d, v[:]))

readvalue(m::AbstractMeas{S,T}, d::Dict) where {S,T} = setvalue(m, getvalue(d, m.tag))
updatethermo(m::AbstractMeas, statevec::AbstractVector{<:Real}, thermo::ThermoModel) = m

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
const VolState{T} = Species{(:V,:T,:P), T, 3} where T
VolState(v::AbstractVector{T}) where T = VolState{T}(v)

@kwdef struct VolumeFlowMeas{S, T, N} <: AbstractSingleMeas{S, T}
    id       :: Symbol
    tag      :: VolState{Union{String,Float64}}
    value    :: VolState{T}
    molarvol :: Species{S, Float64, N}
    stream   :: StreamRef{S, N}
    stdev    :: Float64
end
VolumeFlowMeas{S, T}(x...) where {S,T} = VolumeFlowMeas{S, T, length(S)}(x...)
VolumeFlowMeas{S, T}(;kw...) where {S,T} = VolumeFlowMeas{S, T, length(S)}(;kw...)

function VolumeFlowMeas(measinfo::MeasInfo, streams::Dict{Symbol, <:StreamRef{S}}, thermo::ThermoModel) where S
    if length(measinfo.tags) != 3
        error("Measurement Type: VolumeFlowMeas only supports 3 tags, measurement id '$(measinfo.id)' contains $(length(measinfo.tags))")
    end
    N = length(S)
    (T, P) = (298.15, 101.3e3)
    stream = streams[measinfo.stream]

    thermostate = ThermoState{S, Float64}(
        model=thermo, 
        T=T, 
        P=P, 
        n=Species{S}(ones(N)./N),
        phase=stream.phase
    )

    return VolumeFlowMeas{S, Float64}(
        id       = measinfo.id,
        tag      = VolState{Union{String,Float64}}(measinfo.tags),
        value    = VolState{Float64}(V=0.0, T=T, P=P),
        stdev    = measinfo.stdev[:V],
        stream   = stream,
        molarvol = molar_volumes(thermostate)
    )
end

function readvalue(m::VolumeFlowMeas{S,T}, d::Dict) where {S,T}
    _getvalue(v::String) = T(d[v])
    _getvalue(v::Real)   = T(v)

    return setvalue(m, VolState{T}(_getvalue.(m.tag[:])))
end

function prediction(x::AbstractVector, m::VolumeFlowMeas)
    stream = speciesvec(x, m.stream)
    return dot(m.molarvol[:], stream)
end

function innovation(x::AbstractVector, m::VolumeFlowMeas)
    return m.value[1] - prediction(x, m)
end

function updatethermo(meas::VolumeFlowMeas{S}, statevec::AbstractVector{<:Real}, thermo::ThermoModel) where S
    thermostate = ThermoState(
        model = thermo, 
        T = meas.value[:T], 
        P = meas.value[:P], 
        n = statevec[meas.stream],
        phase = meas.stream.phase
    )
    return @set meas.molarvol = molar_volumes(thermostate)
end

#=============================================================================
Mass flow rates
=============================================================================#
@kwdef struct MassFlowMeas{S, T, N} <: AbstractSingleMeas{S, T}
    id          :: Symbol
    tag         :: String
    value       :: T
    molarmass   :: Species{S, Float64, N}
    stream      :: StreamRef{S, N}
    stdev       :: Float64
end
MassFlowMeas{S, T}(x...) where {S,T}  = MassFlowMeas{S, T, length(S)}(x...)
MassFlowMeas{S, T}(;kw...) where {S,T}  = MassFlowMeas{S, T, length(S)}(;kw...)

function prediction(x::AbstractVector, m::MassFlowMeas)
    stream = speciesvec(x, m.stream)
    return dot(m.molarmass[:], stream[:])
end

function MassFlowMeas(measinfo::MeasInfo, streams::Dict{Symbol, <:StreamRef{S}}, thermo::ThermoModel) where S
    if length(measinfo.tags) != 1
        error("Measurement Type: MassFlowMeas only supports 1 tag, measurement id '$(measinfo.id)' contains $(length(measinfo.tags))")
    end
    
    return MassFlowMeas{S, Float64}(
        id        = measinfo.id,
        tag       = first(values(measinfo.tags)),
        value     = 0.0,
        stdev     = first(values(measinfo.stdev)),
        stream    = streams[measinfo.stream],
        molarmass = molar_weights(thermo)
    )
end

function updatethermo(meas::MassFlowMeas{S}, statevec::AbstractVector{<:Real}, thermo::ThermoModel) where S
    return @set meas.molarmass = molar_weights(thermo)
end

#=============================================================================
Molar Analysis
=============================================================================#
@kwdef struct MoleAnalyzer{S, T, N} <: AbstractMultiMeas{S, T}
    id       :: Symbol
    tag      :: Species{S, String, N}
    value    :: Species{S, T, N}
    stream   :: StreamRef{S, N}
    stdev    :: Species{S, Float64, N}
end
MoleAnalyzer{S, T}(x...) where {S,T} = MoleAnalyzer{S, T, length(S)}(x...)
MoleAnalyzer{S, T}(;kw...) where {S,T} = MoleAnalyzer{S, T, length(S)}(;kw...)

function prediction(x::AbstractVector, m::MoleAnalyzer)
    stream = speciesvec(x, m.stream)
    return stream./sum(stream)
end

function MoleAnalyzer(measinfo::MeasInfo, streams::Dict{Symbol, <:StreamRef{S}}, thermo::ThermoModel) where S
    measid = measinfo.id
    N = length(S)

    if length(measinfo.tags) != N
        error("Measurement Type: MassFlowMeas only supports $(N) tags, measurement id '$(measid)' contains $(length(measinfo.tags))")
    end

    return MoleAnalyzer{S, Float64}(
        id       = measid,
        tag      = Species{S,String}(measinfo.tags),
        value    = zero(Species{S,Float64,N}),
        stdev    = Species{S,Float64,N}(measinfo.stdev),
        stream   = streams[measinfo.stream]
    )
end


#=============================================================================
Mole Balancer
=============================================================================#
@kwdef struct MoleBalance{S, T, N} <: AbstractMultiMeas{S, T}
    id        :: Symbol
    value     :: Species{S, T, N}
    interval  :: Base.RefValue{Float64}
    inlets    :: Vector{StreamRef{S, N}}
    outlets   :: Vector{StreamRef{S, N}}
    reactions :: Vector{ReactionRef{S, N}}
    stdev     :: Species{S, Float64, N}
end
MoleBalance{S, T}(x...) where {S,T} = MoleBalance{S, T, length(S)}(x...)
MoleBalance{S, T}(;kw...) where {S,T} = MoleBalance{S, T, length(S)}(;kw...)

function addinds!(inds::BitVector, m::MoleBalance)
    addinds!(inds, m.inlets)
    addinds!(inds, m.outlets)
    addinds!(inds, m.reactions)
    return inds
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

function MoleBalance(nodeinfo::NodeInfo, streams::Dict{Symbol, <:StreamRef{S}}, interval::Base.RefValue{Float64}) where S
    nodeid = nodeinfo.id
    N = length(S)

    return MoleBalance{S, Float64}(
        id        = nodeid,
        value     = zero(Species{S, Float64, N}),
        interval  = interval,
        stdev     = Species{S}(nodeinfo.stdev),
        inlets    = [streams[id] for id in nodeinfo.inlets],
        outlets   = [streams[id] for id in nodeinfo.outlets],
        reactions = nodeinfo.reactions
    )
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

function populate!(m::MeasCollection, info::MeasInfo, streams::Dict{Symbol, StreamRef})
    type = info.type
    return push!(m[type], type(info, streams))
end

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

function updatethermo!(vmeas::AbstractVector{M}, statevec::AbstractVector{<:Real}, thermo::ThermoModel) where {M <: AbstractMeas}
    updater(m) = updatethermo(m, statevec, thermo)
    vmeas .= updater.(vmeas)
    return vmeas
end

function updatethermo!(c::MeasCollection, statevec::AbstractVector{<:Real}, thermo::ThermoModel)
    for fn in fieldnames(MeasCollection)
        updatethermo!(c[fn], statevec, thermo)
    end
    return c
end

function negloglik(x::AbstractVector, c::MeasCollection)
    return sum(fn-> negloglik(x, c[fn]), fieldnames(MeasCollection))
end
