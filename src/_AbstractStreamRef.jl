include("_ThermoModel.jl")
using Accessors

abstract type AbstractStreamRef{L} end

function stateindex!(indref::Base.RefValue, s::Species{L}) where {L}
    N = length(L)
    start = indref[] + 1
    indref[] = indref[] + N
    return Species{L}(SVector{N}(start:indref[]))
end

function stateindex!(indref::Base.RefValue, r::Integer)
    indref[] = indref[] + 1
    return indref[]
end


#=============================================================================
Construction info for streams
=============================================================================#
@kwdef struct StreamRef{L, N} <: AbstractStreamRef{L}
    id :: Symbol
    index :: Species{L, Int, N}
    phase :: Symbol = :unknown
    refid :: Symbol = :nothing
    scale :: Int = 0
end

function StreamRef{L}(;id, refid=:nothing, phase=:unknown) where L
    N = length(L)
    return StreamRef{L,N}(
        id = id,
        index = Species{L, Int, N}(zero(SVector{N,Int})),
        phase = phase,
        refid = refid,
        scale = 0
    )
end

function stateindex!(indref::Base.RefValue, streaminfo::StreamRef{L}) where {L}
    if streaminfo.refid == :nothing
        return @set streaminfo.index = stateindex!(indref, streaminfo.index)
    else
        return @set streaminfo.scale = stateindex!(indref, streaminfo.scale)
    end
end

function speciesvec(X::AbstractVector, ind::StreamRef)
    species = speciesvec(X, ind.index)
    if ind.refid == :nothing
        return species
    else
        return species.*X[ind.scale]
    end
end

function Base.getindex(X::AbstractVector, ind::StreamRef{L}) where {L}
    return Species{L}(speciesvec(X, ind))
end



"""
fillrefs!(streams::Vector{<:StreamRef})
Once all StremRef values have had their indices set, referenced streams can copy the original streams indices
"""
function fillrefs!(streams::Vector{<:StreamRef})
    streamdict = Dict(stream.id=>stream for stream in streams)

    for (ii, stream) in enumerate(streams)
        if stream.refid != :nothing
            streams[ii] = @set stream.index = streamdict[stream.refid].index 
        end
    end
    return streams
end



#=============================================================================
Chemical reactions
=============================================================================#
@kwdef struct ReactionRef{L,N} <: AbstractStreamRef{L}
    id :: Symbol
    stoich :: Species{L,Float64,N}
    extent :: Int = 0
end

function ReactionRef{L}(;id, stoich) where {L} 
    return ReactionRef{L,length(L)}(
        id = id,
        stoich = stoich,
        extent = 0
    )
end

function stateindex!(indref::Base.RefValue, r::ReactionRef{L}) where {L}
    indref[] = indref[] + 1
    return @set r.extent = indref[]
end

speciesvec(x::AbstractVector, idx::ReactionRef{L}) where L = x[idx.extent] .* speciesvec(idx.stoich)

function Base.getindex(x::AbstractVector{T}, idx::ReactionRef{L}) where {L,T}
    return Species{L}(speciesvec(x, idx))
end


#=============================================================================
Stoichometric relationships between Species vectors and reaction coefficients
=============================================================================#

#Find the extent of reaction based on 
#It is based on which species are consumed (negative values means species is consumed)
#Species that are NOT CONSUMED have an infinite extent (they don't limit the extent of reaction)
function _reagent_extent(rxn_coeff::T1, reagent_moles::T2) where {T1<:Real,T2<:Real} 
    T = promote_type(T1,T2,Float64)
    return ifelse(rxn_coeff<0, reagent_moles/abs(rxn_coeff), T(Inf))
end

"""
stoich_extent(reaction::AbstractVector, input::AbstractVector)

Finds maximum reaction extent based on stoichiometry and the limiting reagent
"""
stoich_extent(stoich::Species{L}, reagents::Species{L}) where L = mapreduce(_reagent_extent, min, stoich[:], reagents[:])
stoich_extent(reaction::ReactionRef{L}, reagents::Species{L}) where L = stoich_extent(reaction.stoich, reagents)

#=============================================================================
Process nodes
=============================================================================#
@kwdef struct NodeRef{L, N}
    id :: Symbol
    stdev     :: Species{L, Float64, N}
    inlets    :: Vector{StreamRef{L, N}}
    outlets   :: Vector{StreamRef{L, N}}
    reactions :: Vector{ReactionRef{L, N}}
end

function stateindex!(indref::Base.RefValue, noderef::NodeRef{L}) where {L}
    for (ii, reaction) in enumerate(noderef.reactions)
        noderef.reactions[ii] = stateindex!(indref, reaction)
    end
    return noderef
end
