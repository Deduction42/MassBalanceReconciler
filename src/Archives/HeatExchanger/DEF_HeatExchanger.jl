abstract type AbstractHeatExchanger{T} <: AbstractAsset{T} end

# ==============================================================================
# General heat exchanger
# ==============================================================================
struct HeatExchangerSpecs <: AbstractSpecs
    mCp :: Float64
    UA  :: Float64
    hot_flow_A_Cf :: Float64
    cold_flow_A_Cf :: Float64
    is_counter_flow::Bool
    is_cross_flow::Bool
end


mutable struct HeatExchanger{T} <: AbstractHeatExchanger{T}
    @inherit_asset_fields(GenericAsset{T})
    LMTD::T
    R_fouling::T
    heating_inlet::ProcessStream{T}
    heating_outlet::ProcessStream{T}
    cooling_inlet::ProcessStream{T}
    cooling_outlet::ProcessStream{T}
    specs::HeatExchangerSpecs
end

function default_value(::Type{asset}, fn::Symbol) where asset <: AbstractHeatExchanger{<:Number}
    T = eltype(asset)
    if fn == :LMTD
        return T(0.1)
    elseif fn == :R_fouling 
        return T(0.0)
    else
        error("Asset type \"$(asset)\" has no default value for field \"$(fn)\"")
    end
end

function default_value(::Type{asset}, fn::Symbol) where asset <: AbstractHeatExchanger{<:Bounds}
    if fn == :LMTD
        return Bounds(0, 1e4)
    elseif fn == :R_fouling 
        return Bounds(0, 1e10)
    else
        error("Asset type \"$(asset)\" has no default value for field \"$(fn)\"")
    end
end



#=
mutable struct HeatExchanger{T,S1,S2,S3,S4} <: AbstractHeatExchanger{T}
    @inherit_asset_fields(GenericAsset{T})
    LMTD::T
    R_fouling::T
    heating_inlet::S1
    heating_outlet::S2
    cooling_inlet::S3
    cooling_outlet::S4
    specs::HeatExchangerSpecs
end
=#