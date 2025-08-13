include("_AbstractStreamRef.jl")

using Accessors
using LinearAlgebra
using FlexUnits, .UnitRegistry
import Base.RefValue


abstract type AbstractMeas{S, T} end
abstract type AbstractSingleMeas{S,T} <: AbstractMeas{S,T} end
abstract type AbstractMultiMeas{S,T} <: AbstractMeas{S,T} end


#==========================================================================================================
Abstract Interface ("value" and "stdev" must be fields)
==========================================================================================================#
eltype(::Type{<:AbstractMeas{S,T}}) where {S,T} = T
meastype(::Type{M}) where M <: AbstractMeas = Base.typename(M).wrapper
meastype(m::AbstractMeas) = meastype(typeof(m))
nan2zero(x::Real) = ifelse(isnan(x), zero(x), x)

#state index collection =======================================================
function addinds!(inds::BitArray, v::AbstractVector{<:Integer})
    inds[v] .= true
    return inds 
end

function addinds!(inds::BitArray, ind::Integer)
    inds[ind] = true
    return inds 
end

function addinds!(inds::BitArray, s::StreamRef) 
    addinds!(inds, s.index[:])
    if hasparents(s)
        addinds!(inds, s.scale)
    end
    return inds
end

function addinds!(inds::BitArray, v::AbstractVector)
    for x in v
        addinds!(inds, x)
    end
    return inds 
end

addinds!(inds::BitArray, r::ReactionRef{L}) where L = addinds!(inds, r.extent)
addinds!(inds::BitArray, m::AbstractMeas) = addinds!(inds, m.stream)


#Tag-reading interface
"""
    getvalue(d::Dict, t::TagInfo)

Uses "t::TagInfo" to look up the tag in "d::Dict", if the tag has no name, return the default value instead
"""
getvalue(d::Dict, t::TagInfo) = hasname(t) ? d[getname(t)] : getvalue(t)
getvalue(d::Dict, v::AbstractVector{TagInfo}) = map(Base.Fix1(getvalue, d), v)
getvalue(d::Dict, v::Species{S, TagInfo}) where S = Species{S}(getvalue(d, v[:]))

"""
    getvalue(d::Dict, t::TagInfo, ind::Integer)

Uses "t::TagInfo" to look up the tag in "d::Dict" at index "ind::Integer", if the tag has no name, return the default value instead
Useful for looking up values when Dict contains arrays of evenly-sampled data
"""
getvalue(d::Dict, t::TagInfo, ind::Integer) = hasname(t) ? d[getname(t)][ind] : getvalue(t)
getvalue(d::Dict, v::AbstractVector{TagInfo}, ind::Integer) = map(t->getvalue(d, t, ind), v)
getvalue(d::Dict, v::Species{S, TagInfo}, ind::Integer) where S = Species{S}(getvalue(d, v[:], ind))

#Retreival functions
"""
    getvalue(m::AbstractMeas{S,T})

Returns the current measured value of "m::AbstractMeas"
"""
getvalue(m::AbstractMeas{S,T}) where {S,T} = m.value
getstdev(m::AbstractMeas{S,T}) where {S,T} = m.stdev
updatethermo!(m::AbstractMeas, statevec::AbstractVector{<:Real}, thermo::ThermoModel) = m
nextindex!(ind::RefValue{<:Integer}, m::AbstractMeas) = nextindices!(ind, valuelength(m))
nextindex!(ind::RefValue{<:Integer}, m::AbstractSingleMeas) = nextindex!(ind)

function setvalues(m::M, v::T) where {T, L, M<:AbstractMeas{L}}
    getfn(fn::Symbol) = (fn == :value) ? v : getproperty(m, fn)
    return basetype(M){L,T}(map(getfn, fieldnames(M))...)
end

"""
    readvalues!(m::AbstractMeas{S,T}, d::Dict, ind::Integer...)

Uses tag information in "m::AbstractMeas{S,T}" to look up values in "d::Dict" and write values to "m"
If an integer is supplied, the dictionary lookup result is indexed (useful if "d" contains arrays)
"""
function readvalues!(m::AbstractMeas{S,T}, d::Dict, ind::Integer...) where {S,T} 
    m.value = getvalue(d, m.tag, ind...)
    return m 
end

function loglik(x::AbstractVector{T}, m::AbstractVector{<:AbstractMeas}) where T <: Real
    RT = promote_type(T,Float64)
    if isempty(m)
        return zero(RT)
    else
        return sum(Base.Fix1(loglik, x), m)
    end
end

function setmeasurement!(v::AbstractVector, ind::RefValue{<:Integer}, m::AbstractMeas)
    v[nextindex!(ind, m)] = getvalue(m)
    return v 
end

function setvariance!(v::AbstractVector, ind::RefValue{<:Integer}, m::AbstractMeas)
    v[nextindex!(ind, m)] = getvariance(m)
    return v
end

function setprediction!(v::AbstractVector, ind::RefValue{<:Integer}, m::AbstractMeas, x::AbstractVector)
    v[nextindex!(ind, m)] = prediction(x, m)
    return v
end

#==========================================================================================================
Univariate measuremnts
==========================================================================================================#
function loglik(x::AbstractVector{T}, m::AbstractSingleMeas) where T <: Real
    return -0.5*abs2(innovation(x, m)/getstdev(m))
end

function innovation(x::AbstractVector, m::AbstractSingleMeas)
    return nan2zero(getvalue(m) - prediction(x, m)) #NaN is considered a missing value, which is ignored
end

function innovation(x::AbstractVector{T}, vm::AbstractVector{AbstractSingleMeas}) where T <: Real
    return map(Base.Fix1(innovation, x), vm)
end

getvariance(m::AbstractSingleMeas) = abs2(getstdev(m))
getnoisecov(m::AbstractSingleMeas) = getvariance(m)
valuelength(m::AbstractSingleMeas) = 1


#==========================================================================================================
Multivariate measuremnts
==========================================================================================================#
getvalue(m::AbstractMultiMeas) = m.value[:]

function loglik(x::AbstractVector{T}, m::AbstractMultiMeas) where T <: Real
    return -0.5*sum(x->x*x, innovation(x, m)./getstdev(m))
end

function innovation(x::AbstractVector, m::AbstractMultiMeas)
    return nan2zero.(getvalue(m) .- prediction(x, m))
end

function innovation(x::AbstractVector{T}, vm::AbstractVector{AbstractMultiMeas{S}}) where {S, T<:Real}
    result = promote_type(T, Float64)[]
    for m in vm
        append!(result, innovation(x,m))
    end
    return result
end

getvariance(m::AbstractMultiMeas) = abs2.(getstdev(m))
getnoisecov(m::AbstractMultiMeas) = Diagonal(getvariance(m))
valuelength(m::AbstractMultiMeas{S}) where {S} = length(S)


#==========================================================================================================
Measurement references (easily find a measurement in a vector with TagInfo)
==========================================================================================================#
@kwdef struct MeasReference
    tag  :: TagInfo
    type :: UnionAll
    ind  :: Int
end

function MeasReference(tag::TagInfo, measvec::AbstractVector{<:AbstractMeas}, measvecs::AbstractVector{<:AbstractMeas}...)
    measref = build_measref(tag, measvec, measvecs...)

    if (measref isa MeasReference)
        return measref 
    else
        throw(ArgumentError("Tag '$(m.flowref.tag)' was not found in any of the provided measurement vectors"))
    end
end


function Base.getindex(v::AbstractVector{<:AbstractMeas}, flowref::MeasReference)
    meas = v[flowref.ind]
    (meas.tag.name == flowref.tag.name) || error("Mismatched Tag: retrieved measurement with tag '$(meas.tag.name)' but was referencing a tag named '$(flowref.tag.name)', measurement vector may have been inappropriately mutated")
    return meas 
end


"""
    build_measref(tag::TagInfo, measvec::AbstractVector{M}) where M <: AbstractMeas

Applies build_measref(tag, ...) to multiple measurement vectors, and returns the first measurement that matches
"""
function build_measref(tag::TagInfo, measvec::AbstractVector{<:AbstractMeas}, measvecs::AbstractVector{<:AbstractMeas}...)
    result = build_measref(tag, measvec)
    return isnothing(result) ? build_measref(tag, measvecs...) : result 
end


"""
    build_measref(tag::TagInfo, measvec::AbstractVector{M}) where M <: AbstractMeas

Builds a measurement reference for a tag, allowing you to easily retrieve that measurement from a vector
This is particularly useful when indexing a MeasCollection, when using this, mesaurement vectors should not 
be mutated other than appending (a check exists to ensure the tag name matches)
"""
function build_measref(tag::TagInfo, measvec::AbstractVector{M}) where M <: AbstractMeas
    ind = findfirst(x->(x.tag.name == tag.name), measvec)
    if isnothing(ind)
        return nothing
    else 
        return MeasReference(tag=tag, type=M.name.wrapper, ind=ind)
    end
end


