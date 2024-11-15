using StaticArrays

#Species inherits all constructors from StaticVector
#Only need to define specific constrcutors as well as supporting Tuple (for varargs)

struct Species{L,T,N} <: StaticVector{N,T}
    data :: SVector{N,T}
    Species{L,T,N}(x::Tuple) where {L,T,N}    = new{L,T,length(L)}(SVector{length(L),T}(x))
    Species{L,T,N}(x::SVector) where {L,T,N}  = new{L,T,length(L)}(SVector{length(L),T}(x))
    Species{L,T}(x::Tuple)   where {L,T}      = new{L,T,length(L)}(SVector{length(L),T}(x))
    Species{L,T}(x::SVector) where {L,T}      = new{L,T,length(L)}(SVector{length(L),T}(x))
    Species{L}(x::SVector{N,T}) where {L,T,N} = new{L,T,length(L)}(SVector{length(L),T}(x))
end

Species{L}(x::Tuple) where {L} = Species{L}(SVector{length(L)}(x))

Species{L,T}(x::NamedTuple) where {L,T} = Species{L,T}(values(x[L]))
Species{L}(x::NamedTuple) where {L}     = Species{L}(values(x[L]))
Species{L,T}(;kwargs...) where {L,T}    = Species{L,T}(kwargs[L])
Species{L}(;kwargs...) where {L}        = Species{L}(kwargs[L])

Base.getindex(x::Species, i::Int) = x.data[i]
function Base.getindex(x::Species{L,T,N}, i::Symbol) where {L,T,N}
    nt = NamedTuple{L}(Base.OneTo(length(L)))
    return x.data[nt[i]]
end

species(::Type{Species{L,T,N}}) where {L,T,N} = L
species(x::Type{Species{L,T,N}}) where {L,T,N} = L
Base.propertynames(x::Type{Species{L,T,N}}) = L