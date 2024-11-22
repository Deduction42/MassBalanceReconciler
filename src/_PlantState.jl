include("_AbstractMeas.jl")

@kwdef struct PlantState{L, F<:Function, N<:Int} 
    statevec :: Vector{Float64}
    statecov :: Matrix{Float64}
    predictor :: F
    measurements :: MeasCollection
end

