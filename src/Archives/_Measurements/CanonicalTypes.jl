# ===============================================================================================
# Generic direct measurements (last resort, slow)
# ===============================================================================================
mutable struct DirectMeasurement{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
    asset_field :: Symbol
end
predict_measurement(m::DirectMeasurement, asset::AbstractAsset) = getproperty(asset, m.asset_field)

# ===============================================================================================
# Basic Stream Measurements
# ===============================================================================================
mutable struct MassFlowMeasurement{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
end
predict_measurement(m::MassFlowMeasurement, asset::ProcessStream) = sum(asset.components)

mutable struct VolFlowMeasurement{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
end
predict_measurement(m::VolFlowMeasurement, asset::ProcessStream) = sum(asset.components .* get_component_volumes(asset))

mutable struct TemperatureMeasurement{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
end
predict_measurement(m::TemperatureMeasurement, asset::AbstractAsset) = asset.temperature

mutable struct PressureMeasurement{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
end
predict_measurement(m::PressureMeasurement, asset::AbstractAsset) = asset.pressure

mutable struct DensityMeasurement{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
end
predict_measurement(m::DensityMeasurement, asset::AbstractAsset)  = get_density(asset)


# ===============================================================================================
# Concentration measurements
# ===============================================================================================
ith_fraction(v::AbstractVector, i::Integer) = v[i]/nonzero(sum(v))

mutable struct MassAnalyzer{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
    component :: Symbol
end

function predict_measurement(m::MassAnalyzer, asset::ProcessStream; idx=COMPONENT_INDEXER)
    return ith_fraction(asset.components, idx[m.component])
end

mutable struct MoleAnalyzer{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
    component :: Symbol
end

function predict_measurement(m::MoleAnalyzer, asset::ProcessStream; idx=COMPONENT_INDEXER)
    components = asset.components ./ CHEM.MolWeight
    return ith_fraction(components, idx[m.component])
end

# ===============================================================================================
# Electrical measurements
# ===============================================================================================
mutable struct Voltmeter{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
end
predict_measurement(m::Voltmeter, asset::AbstractAsset) = asset.voltage

mutable struct Ammeter{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
end
predict_measurement(m::Ammeter, asset::AbstractAsset)   = asset.current

mutable struct Wattmeter{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
end
predict_measurement(m::Wattmeter, asset::AbstractAsset) = asset.voltage * asset.current


# ===============================================================================================
# Mechanical measurements
# ===============================================================================================
mutable struct Tachometer{T} <: AbstractMeasurement{T}
    @inherit_asset_fields(BaseMeasurement{T})
end
predict_measurement(m::Tachometer, asset::AbstractAsset) = asset.angular_velocity