using Clapeyron

include("_Species.jl")

import Clapeyron.molecular_weight

function label2symbol(name::AbstractString)
    return Symbol(replace(lowercase(name), r"\s+"=>"_"))
end

abstract type AbstractThermo{L} end
species(::Type{<:AbstractThermo{L}}) where L = L
species(x::AbstractThermo{L}) where L = L

#=======================================================================================
# Thermodynamic models
=======================================================================================#

@kwdef struct ThermoModel{L, N} <: AbstractThermo{L}
    pure   :: Species{L, PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}, N}
    mixed  :: PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}
    ThermoModel{L,N}(x...) where {L,N} = new{L, length(L)}(x...)
    ThermoModel{L}(x...) where L = new{L, length(L)}(x...)
end

function ThermoModel{L}(thermomap::Dict{Symbol,String})
    species_vec = [thermomap[s] for s in L] 
    models = Species{L}(PR.(species_vec))
    mixed  = PR(species_vec)

    return ThermoModel{L}(models, mixed)
end




#=======================================================================================
# Thermodynamic states
=======================================================================================#
@kwdef struct ThermoState{L,T,N} <: AbstractThermo{L}
    model :: ThermoModel{L,N}
    T :: T
    P :: T 
    n :: Species{L,T,N}
    phase :: Symbol = :unknown
end

molar_weights(model::ThermoModel{L}) where L = Species{L}(molecular_weight.(model.pure))
molar_weights(state::ThermoState) = molar_weights(state.model)

function molar_volumes(state::ThermoState{L}) where L 
    x = state.n[:]./sum(state.n[:])

    mixedvol = volume(state.model.mixed, state.P, state.T, x, phase=state.phase)
    purevol  = volume.(state.model.pure, state.P, state.T, 1.0, phase=state.phase)
    return Species{L}(purevol.*(mixedvol/sum(purevol.*x)))
end




#=======================================================================================
# Test code
=======================================================================================#
const CLAPEYRON_MAP = Dict{Symbol,String}(
    :methane => "methane",
    :ethane => "ethane",
    :propane => "propane",
    :i_butane => "isobutane",
    :n_butane => "butane" ,
    :i_pentane => "isopentane",
    :n_pentane => "pentane",
    :hexane => "hexane",
    :heptane => "heptane",
    :octane  => "octane",
    :nonane => "nonane",
    :decane => "decane",
    :nitrogen => "nitrogen",
    :co2 => "carbon dioxide",
    :oxygen => "oxygen",
    :h2o => "water",
    :co => "carbon monoxide",
    :h2s => "hydrogen sulfide",
    :hydrogen => "hydrogen",
    :helium => "helium",
    :argon => "argon"
)

L = Tuple(collect(keys(CLAPEYRON_MAP)))
model = ThermoModel{L}(CLAPEYRON_MAP)
z = Species{L}(rand(length(L)))
molar_weights(model)
molar_volumes(model, 273.15+30, 101.3*1000, z, phase=:gas)



