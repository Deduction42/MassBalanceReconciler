# ==============================================================================
# General Separation Vessel
# ==============================================================================
abstract type AbstractSeparationVessel <: AbstractAsset end

struct SeparationVesselSpecs <: AbstractSpecs
    diameter  :: Float64
    height    :: Float64
    vol2level :: Array{Float64,2}
end

mutable struct SeparationVessel{T} <: AbstractSeparationVessel
    @inherit_asset_fields(GenericAsset)
    temperature::T
    pressure::T
    
    contents::ProcessStream{T}
    frac_evaporated::Vector{T}

    inlets::Vector{ProcessStream{T}}
    liquid_outlet::ProcessStream{T}
    gas_outlet::ProcessStream{T}
    
    specs::SeparationVesselSpecs
    SeparationVessel{T}() where T = new{T}()
end

# ==============================================================================
# General Separation Tower/Stage (in/out gas, in/out liquid)
# ==============================================================================
abstract type AbstractSeparationStage <: AbstractAsset end
abstract type AbstractSeparationTower <: AbstractSeparationStage end

struct SeparationTowerSpecs <: AbstractSpecs
    diameter  :: Float64
    height    :: Float64
    stages    :: Int64
end

mutable struct SeparationTower{T} <: AbstractSeparationTower
    @inherit_asset_fields(GenericAsset)
    temperature::T
    pressure::T
    
    contents::ProcessStream{T}
    frac_evaporated::Vector{T}

    liquid_inlet::ProcessStream{T}
    gas_inlet::ProcessStream{T}
    liquid_outlet::ProcessStream{T}
    gas_outlet::ProcessStream{T}
    
    specs::SeparationVesselSpecs
    SeparationVessel{T}() where T = new{T}()
end


function predict_contents(vessel::AbstractSeparationVessel)
    contents = copy(vessel.contents.components)
    ΔT = MOD_ENV.ΔT

    for s in vessel.inlets
        contents .= contents .+ ΔT.*s.components
    end 

    contents .= contents .- ΔT.*vessel.liquid_outlet.components
    contents .= contents .- ΔT.*vessel.gas_outlet.components

    return contents
end

function predict_contents(stage::AbstractSeparationStage)
    contents = copy(vessel.contents.components)
    ΔT = MOD_ENV.ΔT

    contents .= contents .+ ΔT.*vessel.liquid_inlet.components
    contents .= contents .+ ΔT.*vessel.gas_inlet.components

    contents .= contents .- ΔT.*vessel.liquid_outlet.components
    contents .= contents .- ΔT.*vessel.gas_outlet.components

    return contents
end

