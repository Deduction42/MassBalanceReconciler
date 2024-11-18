using Clapeyron

include("_Species.jl")

import Clapeyron.molecular_weight

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

function species2symbol(name::AbstractString)
    return Symbol(replace(lowercase(name), r"\s+"=>"_"))
end

@kwdef struct ThermoModel{L, N}
    models   :: Species{L, PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}, N}
    mixed    :: PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}
    ThermoModel{L,N}(x...) where {L,N} = new{L, length(L)}(x...)
    ThermoModel{L}(x...) where L = new{L, length(L)}(x...)
end


function ThermoModel(species::AbstractVector{Symbol})
    L = Tuple(species)

    species_vec = [CLAPEYRON_MAP[s] for s in species] 
    models = Species{L}(PR.(species_vec))
    mixed  = PR(species_vec)

    return ThermoModel{L}(models, mixed)
end

species(::Type{ThermoModel{L,N}}) where {L,N} = L
species(x::ThermoModel{L,N}) where {L,N} = L

function molar_weights(model::ThermoModel{L}) where L
    return Species{L}(molecular_weight.(model.models))
end

function molar_volumes(model::ThermoModel{L}, T::Real, P::Real, z::Species{L}; phase=:unknown) where L 
    zfrac = z./sum(z)
    mixedvol = volume(model.mixed, P, T, zfrac, phase=phase)
    purevol  = volume.(model.models, P, T, zfrac, phase=phase)
    return Species{L}(purevol.*(mixedvol/sum(purevol)))
end

#=======================================================================================
# Test code
=======================================================================================#

model = ThermoModel(collect(keys(CLAPEYRON_MAP)))
L = species(model)
z = Species{L}(rand(length(L)))
molar_weights(model)
molar_volumes(model, 273.15+30, 101.3*1000, z, phase=:gas)



