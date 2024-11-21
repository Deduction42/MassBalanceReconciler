using StaticArrays

#Species inherits all constructors from StaticVector
#Only need to define specific constrcutors as well as supporting Tuple (for varargs)

struct Species{L,T,N} <: StaticVector{N,T}
    data :: SVector{N,T}
    Species{L,T,N}(x::Tuple)  where {L,T,N}  = new{L, T, length(L)}(SVector{length(L),T}(x))
    Species{L,T}(x::Tuple)    where {L,T}    = new{L, T, length(L)}(SVector{length(L),T}(x)) 
    function Species{L}(x::Tuple) where {L}  
        T = Base.promote_typeof(x...)
        return new{L, T, length(L)}(SVector{length(L),T}(x))
    end
end

#Extend key StaticArray interface functions
StaticArrays.length(::Type{<:Species{L}}) where L = length(L)
StaticArrays.construct_type(::Type{Species{L,T}}, x) where {L,T} = StaticArrays.construct_type(Species{L,T,length(L)}, x)
StaticArrays.construct_type(::Type{Species{L}}, x::AbstractArray{T}) where {L,T} = StaticArrays.construct_type(Species{L,T,length(L)}, x)
StaticArrays.construct_type(::Type{Species{L}}, x::Tuple) where {L} = StaticArrays.construct_type(Species{L,Base.promote_typeof(x...),length(L)}, x)

#Add support to NamedTuples and keyword arguments
Species{L,T,N}(x::NamedTuple) where {L,T,N} = Species{L,T}(values(x[L]))
Species{L,T}(x::NamedTuple) where {L,T} = Species{L,T}(values(x[L]))
Species{L}(x::NamedTuple)  where {L}    = Species{L}(values(x[L]))
Species{L,T,N}(;kwargs...)  where {L,T,N}   = Species{L,T}(kwargs[L])
Species{L,T}(;kwargs...)  where {L,T}   = Species{L,T}(kwargs[L])
Species{L}(;kwargs...) where {L}        = Species{L}(kwargs[L])

#Indexing functions
Base.getindex(x::Species, i::Int) = x.data[i]
function Base.getindex(x::Species{L,T,N}, i::Symbol) where {L,T,N}
    nt = NamedTuple{L}(Base.OneTo(length(L)))
    return x.data[nt[i]]
end
Base.getindex(x::Species, i::Colon) = x.data

function Base.get(x::Species{L,T,N}, i::Symbol, d) where {L,T,N}
    nt = NamedTuple{L}(x.data.data)
    return get(nt, i, d)
end
Base.get(x::Species{L,T,N}, i::Nothing, d) where {L,T,N} = d


#Introspection
species(::Type{<:Species{L}}) where L = L
species(x::Species{L}) where L = L
Base.propertynames(x::Species{L}) where L = L

#Populating object based on vector and indexer
function populate(idx::Species{L,<:Integer}, x::AbstractVector{T}) where {L,T} 
    return Species{L,T}(populate_vec(idx, x))
end

function populate_vec(idx::Species{L,<:Integer}, x::AbstractVector{T}) where {L,T} 
    return x[idx.data]
end

#=============================================================================
Chemical reactions
=============================================================================#
struct Reaction{L,T,N}
    extent :: T
    stoich :: Species{L,T,N}
end

Reaction{L,T}(x...) where {L,T} = Reaction{L,T,length(L)}(x...)
Reaction{L}(extent::T, stoich) where {L,T} = Reaction{L,T,length(L)}(extent, stoich)

function populate(idx::Reaction{L,<:Integer}, x::AbstractVector{T}) where {L,T}
    return Species{L}(populate_vec(idx, x))
end

function populate_vec(idx::Reaction{L,<:Integer}, x::AbstractVector{T}) where {L,T}
    return x[idx.extent]*idx.stoich[:]
end