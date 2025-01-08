include("_Species.jl")

using JSON3
abstract type AbstractInfo end

symbolize(x::AbstractString) = Symbol(x)
symbolize(x::AbstractVector{<:AbstractString}) = Symbol.(x)
symbolize(x::AbstractDict{<:AbstractString, T}) where T = Dict(Symbol(k)=>v for (k,v) in pairs(x))
symbolize(x::AbstractDict{<:AbstractString, <:Real}) = Dict(Symbol(k)=>Float64(v) for (k,v) in pairs(x))

symbolize(x::Symbol) = x
symbolize(x::AbstractVector{Symbol}) = Vector(x) 
symbolize(x::AbstractDict{Symbol, T}) where T = Dict(Symbol(k)=>v for (k,v) in pairs(x))
symbolize(x::AbstractDict{Symbol, <:Real}) = Dict(Symbol(k)=>Float64(v) for (k,v) in pairs(x))

symbolize(::Type{T}, x::AbstractDict) where T = Dict{Symbol,T}(Symbol(k)=>convert(T, v) for (k,v) in pairs(x))
