using StaticArrays

struct Species{L,T,N} <: StaticVector{N,T}
    data :: SVector{N,T}
    Species{L,T,N}(x::SVector) where {L,T,N} = new{L,T,N}(x)
end

Species{L,T}(x::SVector{N}) where {L,T,N}  = Species{L,T,length(L)}(x)
Species{L}(x::SVector{N,T}) where {L,T,N}  = Species{L,T,length(L)}(x)

Species{L,T}(x) where {L,T} = Species{L}(SVector{length(L),T}(x))
Species{L}(x) where L = Species{L}(SVector{length(L)}(x))

Species{L,T}(x::AbstractArray) where {L,T}  = Species{L}(SVector{length(L),T}(x))
Species{L}(x::AbstractArray{T}) where {L,T} = Species{L}(SVector{length(L)}(x))

Species{L,T}(x::StaticArray{S}) where {L,T,S} = Species{L}(SVector{length(L),T}(x))
Species{L}(x::StaticArray{S,T}) where {L,T,S} = Species{L}(SVector{length(L)}(x))

Species{L,T}(x::NamedTuple) where {L,T} = Species{L}(SVector{length(L),T}(values(x[L])))
Species{L}(x::NamedTuple) where {L}     = Species{L}(SVector{length(L)}(values(x[L])))

Species{L,T}(;kwargs...) where {L,T} = Species{L,T}(kwargs[L])
Species{L}(;kwargs...) where {L}     = Species{L}(kwargs[L])


Base.getindex(x::Species, i::Int) = x.data[i]
function Base.getindex(x::Species{L,T,N}, i::Symbol) where {L,T,N}
    nt = NamedTuple{L}(Base.OneTo(length(L)))
    return x.data[nt[i]]
end




