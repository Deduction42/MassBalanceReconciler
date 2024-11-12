# ==============================================================================
# Shell and Tube Heat Exchanger
# ==============================================================================
struct HeatExchangerSpecs_ShellAndTube <: AbstractSpecs
    @inherit_spec_fields(HeatExchangerSpecs)
    Passes :: Int64
    NumberOfTubes :: Int64
end

mutable struct HeatExchanger_ShellAndTube{T} <: AbstractHeatExchanger{T}
    @inherit_asset_fields(HeatExchanger{T})
    specs :: HeatExchangerSpecs_ShellAndTube
end

function parent_asset_types(X::HeatExchanger_ShellAndTube{T}) where T
    return [HeatExchanger{T}]
end
