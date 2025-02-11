include("_Species.jl")

using JSON3
using DynamicQuantities
abstract type AbstractInfo end


#Settings for default quantities
const MeasQuantity = Quantity{Float64, SymbolicDimensions{DynamicQuantities.DEFAULT_DIM_BASE_TYPE}}
MeasQuantity(x::String) = parse_units(x)

function parse_units(x::AbstractString)
    ustr = strip(x)
    if isempty(ustr)
        return sym_uparse("1.0") #No units will return a dimensionless 1.0
    else
        return sym_uparse(ustr)
    end
end

function tryparse_units(x::AbstractString) :: Union{String, MeasQuantity}
    try 
        return parse_units(x::AbstractString)
    catch
        return String(x)
    end
end

symbolize(x::AbstractString) = Symbol(x)
symbolize(x::AbstractVector{<:AbstractString}) = Symbol.(x)
symbolize(x::AbstractDict{<:AbstractString, T}) where T = Dict(Symbol(k)=>v for (k,v) in pairs(x))
symbolize(x::AbstractDict{<:AbstractString, <:Real}) = Dict(Symbol(k)=>Float64(v) for (k,v) in pairs(x))

symbolize(x::Symbol) = x
symbolize(x::AbstractVector{Symbol}) = Vector(x) 
symbolize(x::AbstractDict{Symbol, T}) where T = Dict(Symbol(k)=>v for (k,v) in pairs(x))
symbolize(x::AbstractDict{Symbol, <:Real}) = Dict(Symbol(k)=>Float64(v) for (k,v) in pairs(x))

symbolize(::Type{T}, x::AbstractDict) where T = Dict{Symbol,T}(Symbol(k)=>convert(T, v) for (k,v) in pairs(x))
symbolize(::Type{Union{String, MeasQuantity}}, x::AbstractDict{<:AbstractString}) = Dict{Symbol, Union{String, MeasQuantity}}(Symbol(k)=>tryparse_units(v) for (k,v) in pairs(x))
