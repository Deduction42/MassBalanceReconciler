include("_ThermoModel.jl")

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

@kwdef struct ReactionInfo{L, N}
    id :: Symbol
    stoich :: Species{L, Float64, N}
    index  :: Int
end

function populate(rxn::ReactionInfo{L,N}, x::AbstractVector) where {L,N}
    rxnvec = rxn.stoich.data .* x[rxn.index]
    return Species{L, Float64, N}(rxnvec)
end

@kwdef struct NodeInfo{L, T, N}
    id        :: Symbol
    inlets    :: Vector{StreamInfo{L, T, N}}
    outlets   :: Vector{StreamInfo{L, T, N}}
    reactions :: Vector{ReactionInfo{L, T, N}}
end

NodeInfo{L,T}(x...) where {L,T} = NodeInfo{L, T, length(L)}(x...)
NodeInfo{L,T}(;kwargs...) where {L,T} = NodeInfo{L, T, length(L)}(;kwargs...)

function molar_weights(model::ThermoModel{L}, stream::StreamInfo{L,<:Real}) where L
    return molar_weights(model)
end

function molar_volumes(model::ThermoModel{L}, stream::StreamInfo{L,<:Real}) where L
    return molar_volumes(model, stream.t, stream,p, stream.moles, phase=stream.phase)
end
