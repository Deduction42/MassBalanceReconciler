using StaticArrays

##################################################################################################################################
# Chemical species with indexing
##################################################################################################################################

"""
Species{Names::Tuple{Symbol}}

A list of chemical species used to index material streams
    This object behaves like NamedTuple{Names}(Base.OneTo(length(Names)))
"""
struct Species{Names} end

function Base.getindex(::Type{Species{Names}}, i::Symbol) where Names
    nt = NamedTuple{Names}(Base.OneTo(length(Names)))
    return nt[i]
end
Base.length(::Type{Species{Names}}) where Names = length(Names)
species(::Type{Species{Names}}) where Names = Names

##################################################################################################################################
# Main Material Stream API
##################################################################################################################################

abstract type AbstractMaterialStream{Components<:Species, T, N} end

components(stream::AbstractMaterialStream) = stream.comp

Base.getindex(stream::AbstractMaterialStream, i) = components(stream)[i]
function Base.getindex(stream::AbstractMaterialStream{Components, T}, i::Symbol) where {Components, T}
    return stream[Components[i]]
end

Base.eltype(::Type{<:AbstractMaterialStream{Components, T}}) where {Components, T} = T
Base.eltype(stream::T) where T<:AbstractMaterialStream = eltype(T)


##################################################################################################################################
# Specific material streams
##################################################################################################################################
"""
MoleStream{Components, T, N} <: AbstractMaterialStream{Components, T, N}

Material stream where the components are moles
"""
@kwdef struct MoleStream{S<:Species, T, N} <: AbstractMaterialStream{S, T, N}
    comp  :: SVector{N,T}
    phase :: Symbol = :unknown
    id    :: Symbol = :nothing
    MoleStream{S,T,N}(comp, phase, id) where {S,T,N} = new{S,T,N}(comp, phase, id)
    MoleStream{S,T,N}(comp::NamedTuple, phase, id) where {S,T,N} = new{S,T,N}(SVector(values(comp[species(S)])), phase, id)
    MoleStream{S,T}(comp, phase, id) where {S,T} = new{S,T,length(S)}(comp, phase, id)
    MoleStream{S}(comp::AbstractVector{T}, phase, id) where {S,T} = new{S,T,length(S)}(comp, phase, id)
end

MoleStream{S,T}(;kwargs...) where {S,T} = MoleStream{S,T,length(S)}(;kwargs...)


"""
MassStream{Components, T, N} <: AbstractMaterialStream{Components, T, N}

Material stream where the components are mass
"""
@kwdef struct MassStream{S<:Species, T, N} <: AbstractMaterialStream{S, T, N}
    comp  :: SVector{N,T}
    phase :: Symbol = :unknown
    id    :: Symbol = :nothing
    MassStream{S,T,N}(comp, phase, id) where {S,T,N} = new{S,T,N}(comp, phase, id)
    MassStream{S,T,N}(comp::NamedTuple, phase, id) where {S,T,N} = new{S,T,N}(SVector(values(comp[species(S)])), phase, id)
    MassStream{S,T}(comp, phase, id) where {S,T} = new{S,T,length(S)}(comp, phase, id)
    MassStream{S}(comp::AbstractVector{T}, phase, id) where {S,T} = new{S,T,length(S)}(comp, phase, id)
end

MassStream{S,T}(;kwargs...) where {S,T} = MassStream{S,T,length(S)}(;kwargs...)

##################################################################################################################################
# Material Stream Conversions
##################################################################################################################################
function MassStream(stream::MoleStream{S, T, N}, molwt::AbstractMaterialStream{S,<:Real}) where {S, T, N}
    return MassStream{S, T, N}(
        id = stream.id,
        comp = components(stream).*components(molwt),
        phase = stream.phase
    )
end

function MoleStream(stream::MassStream{S, T, N}, molwt::AbstractMaterialStream{S,<:Real}) where {S, T, N}
    return MoleStream{S, T, N}(
        id = stream.id,
        comp = components(stream)./components(molwt),
        phase = stream.phase
    )
end