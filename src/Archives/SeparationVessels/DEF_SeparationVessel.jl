abstract type AbstractSeparatorVLE{T} <: AbstractAsset{T} end

# ==============================================================================
# General separation vessel
# ==============================================================================
struct SeparatorVLE_Specs <: AbstractSpecs
    Diameter :: Vector{Float64}
    Length   :: Vector{Float64}
    IsVertical :: Bool
end

mutable struct SeparatorVLE{T} <: AbstractSeparatorVLE{T}
    @inherit_asset_fields(GenericAsset{T})
    accumulation    :: T
    heat_input      :: T
    vapor_frac      :: Vector{T}
    contents        :: ProcessStream{T}
    inlet           :: ProcessStream{T}
    gas_outlet      :: ProcessStream{T}
    liquid_outlet   :: ProcessStream{T}
    specs           :: SeparatorVLE_Specs
end

function default_value(::Type{asset}, fn::Symbol) where asset <: AbstractSeparatorVLE{<:Number}
    T = eltype(asset)
    if fn == :accumulation
        return T(0.0)
    elseif fn == :heat_input
        return T(0.1)
    elseif fn == :vapor_frac
        return T.(fill(0.5, length(CHEM.Names)))
    else
        error("Asset type \"$(asset)\" has no default value for field \"$(fn)\"")
    end
end

function default_value(::Type{asset}, fn::Symbol) where asset <: AbstractSeparatorVLE{<:Bounds}
    if fn == :accumulation
        return Bounds(0.0, 1e12)
    elseif fn == :heat_input
        return Bounds(0.0, 1e12)
    elseif fn == :vapor_frac
        return fill(Bounds(0.0, 1.0), length(CHEM.Names))
    else
        error("Asset type \"$(asset)\" has no default value for field \"$(fn)\"")
    end
end