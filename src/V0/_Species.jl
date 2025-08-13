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
speciesvec(x::Species) = x.data

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

#Support for constructing with dictionaries
Species{L,T,N}(d::Dict{Symbol}) where {L,T,N}    = Species{L,T}([d[k] for k in L])
Species{L,T,N}(d::Dict{Symbol}, x) where {L,T,N} = Species{L,T}([get(d, k, x) for k in L])
Species{L,T}(d::Dict{Symbol}) where {L,T}        = Species{L,T}([d[k] for k in L])
Species{L,T}(d::Dict{Symbol}, x) where {L,T}     = Species{L,T}([get(d, k, x) for k in L])
Species{L}(d::Dict{Symbol,T}) where {L,T}        = Species{L,T}([d[k] for k in L])
Species{L}(d::Dict{Symbol,T}, x) where {L,T}     = Species{L,T}([get(d, k, T(x)) for k in L])


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


#Conversion to dictionary
Base.Dict(x::Species{L}) where L = Dict(L .=> speciesvec(x))


#=============================================================================
Interactions with abstract vectors
=============================================================================#
function Base.getindex(x::AbstractVector{T}, idx::Species{L,<:Integer}) where {L,T}
    return Species{L}(x[idx.data])
end

function Base.setindex!(x::AbstractVector{T}, idx::Species{L,<:Integer}, vals::Species{L}) where {L,T}
    x[idx.data] = vals.data
    return x 
end

speciesvec(x::AbstractVector, idx::Species)  = x[idx.data]

#=============================================================================
Default methods for total and specific aggregation (only works for same type)
=============================================================================#
totals(::Type{<:Species{L}}, s::Species{L}) where L = return s
molaravgs(::Type{<:Species{L}}, s::Species{L}, fracs=nothing) where L = return s

#=============================================================================
More infomrative errors prompting users to define their own aggregation rules
=============================================================================#
function totals(::Type{<:Species{L}}, s::Species{S}) where {L,S}
    return throw(ArgumentError("""
    No defined way to use 'totals' (totalizing over applicable components) to convert species 
      'S=$(S)' 
    to 
      'L=$(L)'
    It is up to the user to define such conversion rules by defining:
      totals(::Type{<:Species{L}}, s::Species{S})
    """))
end

function molaravgs(::Type{<:Species{L}}, s::Species{S}, fracs=nothing) where {L,S}
    return throw(ArgumentError("""
    No defined way to use 'molaravgs' (molar averages over applicable components) to convert species
      'S=$(S)' 
    to 
      'L=$(L)'
    It is up to the user to define such conversion rules by defining:
      molaravgs(::Type{<:Species{L}}, s::Species{S}, moles::Species{S})
    """))
end
