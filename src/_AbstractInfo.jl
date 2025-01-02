include("_Species.jl")

abstract type AbstractInfo end

symbolize(x::AbstractString) = Symbol(x)
symbolize(x::AbstractVector{<:AbstractString}) = Symbol.(x)
symbolize(x::AbstractDict{<:AbstractString, T}) where T = Dict{Symbol,T}(Symbol(k)=>v for (k,v) in pairs(x))
symbolize(x::Symbol) = x
symbolize(x::AbstractVector{Symbol}) = x 
symbolize(x::AbstractDict{Symbol, T}) where T = x