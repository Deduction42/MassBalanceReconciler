
abstract type Measurement{S, T} end


@kwdef struct VolumeFlowMeas{S, T, N} <: Measurement{S, T}
    value    :: T
    molarvol :: Species{S, Float64, N}
    flowind  :: Species{S, Int, N}
end
VolumeFlowMeas{S, T}(x...) = VolumeFlowMeas{S, T, length(S)}(x...)

function innovation(x::AbstractVector, m::VolumeFlowMeas)
    stream = populate(m.flowind, x)
    return m.value - dot(species(m.molarvol), species(stream))
end

stateindex(m::VolumeFlowMeas) = collect(m.flowind)


@kwdef struct MassFlowMeas{S, T, N} <: Measurement{T}
    value       :: T
    molarmass   :: Species{S, Float64, N}
    flowind     :: Species{S, Int, N}
end
MassFlowMeas{S, T}(x...) = MassFlowMeas{S, T, length(S)}(x...)