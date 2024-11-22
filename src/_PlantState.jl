include("_AbstractMeas.jl")

@kwdef struct PlantState{L, F, N<:Int} 
    statevec     :: Vector{Float64}
    statecov     :: Matrix{Float64}
    dpredictor   :: Matrix{Float64}
    measurements :: MeasCollection{L, Float64, N}
    streams      :: Vector{StreamInfo{L,N}}
    nodes        :: Vector{NodeInfo{L,N}}
end

