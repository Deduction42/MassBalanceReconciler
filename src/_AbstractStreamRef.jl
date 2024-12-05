include("_ThermoModel.jl")
using Accessors

abstract type AbstractStreamRef{L} end

#=============================================================================
Construction info for streams
=============================================================================#
@kwdef struct StreamRef{L, N} <: AbstractStreamRef{L}
    id :: Symbol
    massflow :: Float64
    index :: Species{L, Int, N}
    phase :: Symbol = :unknown
    refid :: Symbol = :nothing
    scale :: Int = 0
end

function StreamRef{L}(;id, massflow, refid=:nothing, phase=:unknown) where L
    N = length(L)
    return StreamRef{L,N}(
        id = id,
        index = Species{L, Int, N}(zero(SVector{N,Int})),
        massflow = massflow,
        phase = phase
    )
end

function stateindex!(indref::Base.RefValue, streaminfo::StreamRef{L}) where {L}
    if streaminfo.refid == :nothing
        return @set streaminfo.index = stateindex!(indref, streaminfo.index)
    else
        return @set streaminfo.scale = stateindex!(indref, streaminfo.scale)
    end
end

function Base.getindex(X::AbstractVector, ind::StreamRef{L}) where {L}
    return Species{L}(speciesvec(X, ind))
end

function speciesvec(X::AbstractVector, ind::StreamRef)
    species = speciesvec(X, ind.index)
    if ind.refid == :nothing
        return species
    else
        return species.*X[ind.scale]
    end
end






#=============================================================================
Chemical reactions
=============================================================================#
struct Reaction{L,T,N}
    extent :: T
    stoich :: Species{L,T,N}
end
speciesvec(r::Reaction) = r.extent .* r.stoich[:]

Reaction{L,T}(x...) where {L,T} = Reaction{L,T,length(L)}(x...)
Reaction{L}(extent::T, stoich) where {L,T} = Reaction{L,T,length(L)}(extent, stoich)


function Base.getindex(x::AbstractVector{T}, idx::Reaction{L,<:Integer}) where {L,T}
    return Reaction{L}(x[idx.extent], idx.stoich)
end

function Base.setindex!(x::AbstractVector{T}, idx::Reaction{L,<:Integer}, val::Reaction{L}) where {L,T}
    x[idx.extent] = val.extent
    return x
end

speciesvec(x::AbstractVector, idx::Reaction) = x[idx.extent] .* idx.stoich.data

#=============================================================================
Stoichometric relationships between Species vectors and reaction coefficients
=============================================================================#

#Find the extent of reaction based on 
#It is based on which species are consumed (negative values means species is consumed)
#Species that are NOT CONSUMED have an infinite extent (they don't limit the extent of reaction)
function reagent_extent(rxn_coeff::T1, reagent_moles::T2) where {T1<:Real,T2<:Real} 
    T = promote_type(T1,T2,Float64)
    return ifelse(rxn_coeff<0, reagent_moles/abs(rxn_coeff), T(Inf))
end

"""
stoich_extent(reaction::AbstractVector, input::AbstractVector)

Finds maximum reaction extent based on stoichiometry and the limiting reagent
"""
stoich_extent(reaction::AbstractVector, input::AbstractVector) = mapreduce(reagent_extent, min, reaction, input)

