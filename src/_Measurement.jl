
abstract type Measurement{T} end

@kwdef struct VolumeFlowMeas{S, N}
    value    :: Float64
    molarvol :: MoleStream{S, Float64, N}
    stream   :: MoleStream{S, Int64, N}
end

function innovation(x::AbstractVector, m::VolumeFlowMeas)
    stream = populate(m.stream, x)
    return m.value - dot(species(m.molarvol), species(stream))
end