include("_StreamInfo.jl")

abstract type AbstractMeas{S, T} end

stateindex(m::AbstractMeas) = collect(m.flowind)
standarderr(x::AbstractVector, m::AbstractMeas) = innovation(x,m)/m.stdev

#=============================================================================
Volumetric flow rates
=============================================================================#
@kwdef struct VolumeFlowMeas{S, T, N} <: AbstractMeas{S, T}
    value    :: T
    molarvol :: Species{S, Float64, N}
    flowind  :: Species{S, Int, N}
    stdev    :: Float64
end
VolumeFlowMeas{S, T}(x...) where {S,T} = VolumeFlowMeas{S, T, length(S)}(x...)

function innovation(x::AbstractVector, m::VolumeFlowMeas)
    stream = populate(m.flowind, x)
    return m.value - dot(m.molarvol, stream)
end

#=============================================================================
Mass flow rates
=============================================================================#
@kwdef struct MassFlowMeas{S, T, N} <: AbstractMeas{S, T}
    value       :: T
    molarmass   :: Species{S, Float64, N}
    flowind     :: Species{S, Int, N}
    stdev       :: Float64
end
MassFlowMeas{S, T}(x...) where {S,T}  = MassFlowMeas{S, T, length(S)}(x...)

function innovation(x::AbstractVector, m::MassFlowMeas)
    stream = populate(m.flowind, x)
    return m.value - dot(m.molamass, stream)
end

#=============================================================================
Molar Analysis
=============================================================================#
@kwdef struct MoleAnalyzer{S, T, N} <: AbstractMeas{S, T}
    values  :: Species{S,T,N}
    flowind :: Species{S,T,N}
    stdev   :: Float64
end
MoleAnalyzer{S, T}(x...) where {S,T} = MoleAnalyzer{S, T, length(S)}(x...)

function innovation(x::AbstractVector, m::MoleAnalyzer)
    stream = populate(m.flowind, x)
    return m.values .- stream./sum(stream)
end

