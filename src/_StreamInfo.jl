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
    index :: Reaction{L, Int, N}
end

function populate(rxn::ReactionInfo{L,N}, x::AbstractVector) where {L,N}
    rxnvec = rxn.stoich.data .* x[rxn.index]
    return Species{L, Float64, N}(rxnvec)
end

@kwdef struct BalanceInfo{L, T, N}
    id        :: Symbol
    inlets    :: Vector{StreamInfo{L, T, N}}
    outlets   :: Vector{StreamInfo{L, T, N}}
    reactions :: Vector{ReactionInfo{L, N}}
end

BalanceInfo{L,T}(x...) where {L,T} = BalanceInfo{L, T, length(L)}(x...)
BalanceInfo{L,T}(;kwargs...) where {L,T} = BalanceInfo{L, T, length(L)}(;kwargs...)

function molar_weights(model::ThermoModel{L}, stream::StreamInfo{L,<:Real}) where L
    return molar_weights(model)
end

function molar_volumes(model::ThermoModel{L}, stream::StreamInfo{L,<:Real}) where L
    return molar_volumes(model, stream.t, stream,p, stream.moles, phase=stream.phase)
end
