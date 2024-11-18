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
Species{L,T}(x::NamedTuple) where {L,T} = Species{L,T}(values(x[L]))
Species{L}(x::NamedTuple)  where {L}    = Species{L}(values(x[L]))
Species{L,T}(;kwargs...)  where {L,T}   = Species{L,T,N}(kwargs[L])
Species{L}(;kwargs...) where {L}        = Species{L}(kwargs[L])

#Indexing functions
Base.getindex(x::Species, i::Int) = x.data[i]
function Base.getindex(x::Species{L,T,N}, i::Symbol) where {L,T,N}
    nt = NamedTuple{L}(Base.OneTo(length(L)))
    return x.data[nt[i]]
end

#Introspection
species(::Type{Species{L,T,N}}) where {L,T,N} = L
species(x::Type{Species{L,T,N}}) where {L,T,N} = L
Base.propertynames(x::Species{L,T,N}) where {L,T,N} = L