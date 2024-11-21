include("_ThermoModel.jl")

#=============================================================================
Construction info for streams
=============================================================================#
@kwdef struct StreamInfo{L, T, N}
    id :: Symbol
    t :: T
    p :: T
    moles :: Species{L, T, N}
    index :: Species{L, Int, N}
    phase :: Symbol = :unknown
end

StreamInfo{L,T}(x...) where {L,T} = StreamInfo{L, T, length(L)}(x...)
StreamInfo{L,T}(;kwargs...) where {L,T} = StreamInfo{L, T, length(L)}(;kwargs...)

function molar_weights(model::ThermoModel{L}, stream::StreamInfo{L,<:Real}) where L
    return molar_weights(model)
end

function molar_volumes(model::ThermoModel{L}, stream::StreamInfo{L,<:Real}) where L
    return molar_volumes(model, stream.t, stream,p, stream.moles, phase=stream.phase)
end



@kwdef struct NodeInfo{L,N}
    id        :: Symbol
    inlets    :: Vector{Symbol}
    outlets   :: Vector{Symbol}
    reactions :: Vector{Reaction{L, Int, N}}
end

NodeInfo{L}(x...) where {L} = NodeInfo{L,length(L)}(x...)
NodeInfo{L}(;kwargs...) where {L} = NodeInfo{L, length(L)}(;kwargs...)


