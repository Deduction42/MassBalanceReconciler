using Accessors
include("_StreamInfo.jl")


#=============================================================================
Abstract Interface ("value" and "stdev" must be fields)
=============================================================================#
abstract type AbstractMeas{S, T} end

eltype(::Type{<:AbstractMeas{S,T}}) where {S,T} = T

stateindex(m::AbstractMeas) = collect(m.flowind)
standarderr(x::AbstractVector, m::AbstractMeas) = innovation(x,m)/m.stdev

meastype(::Type{M}) where M <: AbstractMeas = Base.typename(M).wrapper
meastype(m::AbstractMeas) = meastype(typeof(m))

changevalue(m::AbstractMeas, v) = @set m.value = v
getvalue(d::Dict{T}, v::T) where T = d[v]
getvalue(d::Dict{T}, v::AbstractVector{T}) where T  = map(Base.Fix1(get, d), v)
getvalue(d::Dict{T}, v::Species{S,T}) where {S,T}   = Species{S}(getvalue(d, v[:]))

readvalue(d::Dict{T}, m::AbstractMeas{S,T}) where {S,T} = changevalue(m, getvalue(d, m.value))

function negloglik(x::AbstractVector{T}, m::AbstractVector{<:AbstractMeas}) where T <: Real
    return sum(Base.Fix1(negloglik, x), m, init=zero(promotetype(T,Float64)))
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

function innovation(x::AbstractVector{T}, vm::AbstractVector{AbstractSingleMeas}) where T <: Real
    result = promote_type(T, Float64)[]
    for m in vm
        append!(result, innovation(x,m))
    end
    return result
end

#=============================================================================
Volumetric flow rates
=============================================================================#
@kwdef struct VolumeFlowMeas{S, T, N} <: AbstractSingleMeas{S, T}
    value    :: T
    molarvol :: Species{S, Float64, N}
    flowind  :: Species{S, Int, N}
    stdev    :: Float64
end
VolumeFlowMeas{S, T}(x...) where {S,T} = VolumeFlowMeas{S, T, length(S)}(x...)

function prediction(x::AbstractVector, m::VolumeFlowMeas)
    stream = populate_vec(m.flowind, x)
    return dot(m.molarvol[:], stream)
end

#=============================================================================
Mass flow rates
=============================================================================#
@kwdef struct MassFlowMeas{S, T, N} <: AbstractSingleMeas{S, T}
    value       :: T
    molarmass   :: Species{S, Float64, N}
    flowind     :: Species{S, Int, N}
    stdev       :: Float64
end
MassFlowMeas{S, T}(x...) where {S,T}  = MassFlowMeas{S, T, length(S)}(x...)

function prediction(x::AbstractVector, m::MassFlowMeas)
    stream = populate_vec(m.flowind, x)
    return dot(m.molarmass[:], stream)
end

#=============================================================================
Molar Analysis
=============================================================================#
@kwdef struct MoleAnalyzer{S, T, N} <: AbstractMeas{S, T}
    value   :: Species{S, T, N}
    flowind :: Species{S, Int, N}
    stdev   :: Species{S, Float64, N}
end
MoleAnalyzer{S, T}(x...) where {S,T} = MoleAnalyzer{S, T, length(S)}(x...)

function prediction(x::AbstractVector, m::MoleAnalyzer)
    stream = populate_vec(m.flowind, x)
    return stream./sum(stream)
end


#=============================================================================
Mole Balancer
=============================================================================#
@kwdef struct MoleBalance{S, T, N} <: AbstractMeas{S, T}
    value     :: Species{S, T, N}
    inlets    :: Vector{Species{S, Int, N}}
    outlets   :: Vector{Species{S, Int, N}}
    reactions :: Vector{Reaction{S, Int, N}}
    stdev     :: Species{S, Float64, N}
end
MoleBalance{S, T}(x...) where {S,T} = MoleBalance{S, T, length(S)}(x...)

function prediction(x::AbstractVector{T}, m::MoleBalance{S, <:Integer, N}) where {S,T,N}
    balinit = zero(SVector{N, promote_type(T,Float64)})

    balance = (
          sum(Base.Fix2(populate_vec, x), m.inlets, init=balinit)
        - sum(Base.Fix2(populate_vec, x), m.outlets, init=balinit)
        + sum(Base.Fix2(populate_vec, x), m.reactions, init=balinit)
    )

    return balance
end

#=============================================================================
Collection of all measurements
=============================================================================#
@kwdef struct MeasCollection
    VolumeFlowMeas  :: Vector{VolumeFlowMeas}   = VolumeFlowMeas[]
    MassFlowMeas    :: Vector{MassFlowMeas}     = MassFlowMeas[]
    MoleAnalyzer    :: Vector{MoleAnalyzer}     = MoleAnalyzer[]
    MoleBalance     :: Vector{MoleBalance}      = MoleBalance[]
end

Base.getindex(m::MeasCollection, k::Symbol)  = getproperty(m, k)
Base.getindex(m::MeasCollection, k::Integer) = getproperty(m, fieldnames(MeasCollection)[k])
Base.getindex(m::MeasCollection, k::Colon)   = map(Base.Fix1(getproperty, m), fieldnames(MeasCollection))
Base.firstindex(m::MeasCollection) = 1
Base.lastindex(m::MeasCollection) = length(fieldnames(MeasCollection))
Base.getindex(m::MeasCollection, ::Type{T}) where T = getproperty(m, Symbol(T))
Base.getindex(m::MeasCollection, k::AbstractVector) = map(Base.Fix1(getindex, m), k)
Base.getindex(m::MeasCollection, k::Tuple) = map(Base.Fix1(getindex, m), k)
