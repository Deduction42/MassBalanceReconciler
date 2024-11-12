# ==============================================================================
# Shell and Tube Heat Exchanger
# ==============================================================================
function parent_asset_types(X::HeatExchanger_ShellAndTube{T}) where T
    return [HeatExchanger{T}]
end

function predict_state!(hx1::HX, hx0::HX, ::Type{HX_Type}) where {HX<:AbstractHeatExchanger, HX_Type<:HeatExchanger_ShellAndTube}
    predict_state!(hx1, hx0, HeatExchanger)
end
