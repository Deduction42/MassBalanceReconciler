include("_ThermoModel.jl")
using Accessors

#=============================================================================
Construction info for streams
=============================================================================#
@kwdef struct StreamInfo <: AbstractInfo
    id        :: Symbol
    massflow  :: Float64
    molefracs :: Union{Symbol, Dict{Symbol, Float64}}
    phase     :: Symbol = :unknown
end
hasparent(info::StreamInfo) = info.molefracs isa Symbol

function StreamInfo(d::AbstractDict{Symbol})
    molefracs = d[:molefracs]

    return StreamInfo(
        id = Symbol(d[:id]),
        massflow  = d[:massflow],
        molefracs = (molefracs isa AbstractDict) ? symbolize(Float64, molefracs) : symbolize(molefracs),
        phase = Symbol(get(d, :phase, :unknown))
    )
end

StreamInfo(d::AbstractDict{<:AbstractString}) = StreamInfo(symbolize(d))

#=============================================================================
Abstract stream interface
=============================================================================#
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
hasparent(streamref::StreamRef) = (streamref.refid != :nothing)

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

function StreamRef{L}(info::StreamInfo) where L
    refid = hasparent(info) ? info.molefracs : :nothing
    return StreamRef{L}(id=info.id, refid=refid, phase=info.phase)
end

function stateindex!(indref::Base.RefValue, streamref::StreamRef{L}) where {L}
    if hasparent(streamref)
        return @set streamref.scale = stateindex!(indref, streamref.scale)
    else
        return @set streamref.index = stateindex!(indref, streamref.index)
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
stateindex!(indref::Base.RefValue, streams::Vector{<:StreamRef})

Fills out state index information for `streams`, retuns the reference of the indexer (holding the last state index)
"""
function stateindex!(indref::Base.RefValue, streams::Vector{<:StreamRef})
    #Fill out state indexing information
    for (ii, stream) in enumerate(streams)
        streams[ii] = stateindex!(indref, stream)
    end

    #Fill all compositions that are missing (because they refer to another stream)
    streamdict = Dict(stream.id=>stream for stream in streams)
    for k in keys(streamdict)
        _fillcomposition!(streamdict, k)
    end

    #Fill the stream list
    for (ii, stream) in enumerate(streams)
        if hasparent(stream)
            streams[ii] = streamdict[stream.id]
        end
    end

    return indref
end

#Fill the composition of the stream indexed by "id" (fills in the parent if it's not filled in yet)
function _fillcomposition!(streamdict::Dict{Symbol, <:StreamRef}, id::Symbol)
    stream = streamdict[id]

    #If the first index is non-zero, this stream composition is already filled out so just return it
    if !iszero(first(stream.index)) 
        return stream
    end

    #Ensure that the parent stream is already filled out
    parent = _fillcomposition!(streamdict, stream.refid)

    #Return the stream with the new index
    newstream = @set stream.index = parent.index 
    streamdict[id] = newstream
    
    return newstream
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

function ReactionRef{L}(id::Symbol, reactinfo::Dict{Symbol,Float64}) where L
    stoichvec = [reactinfo[l] for l in L]
    return ReactionRef{L}(id=id, stoich=Species{L}(stoichvec))
end


function stateindex!(indref::Base.RefValue, r::ReactionRef{L}) where {L}
    indref[] = indref[] + 1
    return @set r.extent = indref[]
end

speciesvec(x::AbstractVector, idx::ReactionRef{L}) where L = x[idx.extent] .* speciesvec(idx.stoich)

function Base.getindex(x::AbstractVector{T}, idx::ReactionRef{L}) where {L,T}
    return Species{L}(speciesvec(x, idx))
end

"""
stateindex!(indref::Base.RefValue, reactions::Vector{<:ReactionRef})

Fills out state index information for `reactions`, retuns the reference of the indexer (holding the last state index)
"""
function stateindex!(indref::Base.RefValue, reactions::Vector{<:ReactionRef})
    for (ii, reaction) in enumerate(reactions)
        reactions[ii] = stateindex!(indref, reaction)
    end
    return indref
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
Construction info for nodes
=============================================================================#
@kwdef struct NodeInfo <: AbstractInfo
    id        :: Symbol
    stdev     :: Dict{Symbol, Float64}
    inlets    :: Vector{Symbol}
    outlets   :: Vector{Symbol}
    reactions :: Vector{Dict{Symbol, Float64}} = Dict{Symbol, Float64}[]
end

function NodeInfo(d::AbstractDict{Symbol})
    return NodeInfo(
        id = Symbol(d[:id]),
        stdev  = symbolize(Float64, d[:stdev]),
        inlets = symbolize(d[:inlets]),
        outlets = symbolize(d[:outlets]),
        reactions = symbolize.(Float64, d[:reactions])
    )
end

NodeInfo(d::AbstractDict{<:AbstractString}) = NodeInfo(symbolize(d))

function add_reaction!(nodeinfo::NodeInfo, stoich::Species)
    push!(nodeinfo.reactions, ReactionRef{L,N}(0, stoich))
end


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

function NodeRef{L}(info::NodeInfo, streamdict::Dict{Symbol, <:StreamRef}) where L
    return NodeRef{L, length(L)}(
        id        = info.id,
        stdev     = Species{L}(info.stdev),
        inlets    = [streamdict[id] for id in info.inlets],
        outlets   = [streamdict[id] for id in info.outlets],
        reactions = [ReactionRef{L}(id=info.id, stoich) for stoich in info.reactions]
    )
end



function stateindex!(indref::Base.RefValue, noderef::NodeRef{L}) where {L}
    return stateindex!(indref, noderef.reactions)
end

function stateindex!(indref::Base.RefValue, noderefs::Vector{<:NodeRef})
    for noderef in noderefs
        stateindex!(indref, noderef)
    end
    return indref
end

#=============================================================================
Construction info for simple stream relationshps, useful for predictions
=============================================================================#
@kwdef struct StreamRelationship
    id     :: Symbol
    parent :: Symbol
    factor :: Float64
    timeconst :: Float64
end

function StreamRelationship(d::AbstractDict{Symbol})
    return StreamRelationship(
        id = Symbol(d[:id]),
        parent = Symbol(d[:parent]),
        factor = d[:factor],
        timeconst = d[:timeconst]
    )
end

function state_transition(Nx::Int, relationships::Vector{StreamRelationship}, streamdict::Dict{Symbol, <:StreamRef})
    #Build the predictor based on relationships, and the noise intensity based off initial state covariance
    transmat = zeros(Nx, Nx)
    for relationship in relationships
        streamind = streamdict[relationship.id].index.data 
        parentind = streamdict[relationship.parent].index.data 
        
        for (s,p) in zip(streamind, parentind)
            transmat[s,s] = -1/relationship.timeconst
            transmat[s,p] = relationship.factor/relationship.timeconst
        end
    end

    return transmat
end