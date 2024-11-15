
abstract type Measurement{T} end

@kwdef struct VolumeFlowMeas{S, T, N} <: Measurement{T}
    value    :: Float64
    molarvol :: Species{S, Float64, N}
    flowind  :: Species{S, Int, N}
end
VolumeFlowMeas{S, T}(x...) = VolumeFlowMeas{S, T, length(S)}(x...)

function innovation(x::AbstractVector, m::VolumeFlowMeas)
    stream = populate(m.flowind, x)
    return m.value - dot(species(m.molarvol), species(stream))
end

function stateindex(m::VolumeFlowMeas)
    return collect(m.flowind)
end