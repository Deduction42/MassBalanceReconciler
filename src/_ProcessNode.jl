include("_MaterialStream.jl")

abstract type AbstractProcessNode{S<:Species, T} end

@kwdef struct ProcessNode{S, T, N}
    inlets  :: Vector{MoleStream{S, T, N}}
    outlets :: Vector{MoleStream{S, T, N}}
end
