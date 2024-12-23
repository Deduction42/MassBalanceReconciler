using Clapeyron

include("_Species.jl")

import Clapeyron.molecular_weight

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

function ThermoModel{L}(thermomap::Dict{Symbol,String}) where {L}
    species_vec = [thermomap[s] for s in L] 
    models = Species{L}(PR.(species_vec))
    mixed  = PR(species_vec)

    return ThermoModel{L}(models, mixed)
end




#=======================================================================================
# Thermodynamic states
=======================================================================================#
@kwdef struct ThermoState{L,ET,N} <: AbstractThermo{L}
    model :: ThermoModel{L,N}
    T :: ET
    P :: ET 
    n :: Species{L,ET,N}
    phase :: Symbol = :unknown
end

ThermoState{L,T}(x...) where {L,T} = ThermoState{L,T,length(L)}(x...)
ThermoState{L,T}(;kwargs...) where {L,T} = ThermoState{L,T,length(L)}(kwargs[fieldnames(ThermoState)]...)

molar_weights(model::ThermoModel{L}) where L = Species{L}(molecular_weight.(model.pure))
molar_weights(state::ThermoState) = molar_weights(state.model)
molar_weights(::Type{<:Species{L}}, state::ThermoState) where L = molaravgs(Species{L}, molar_weights(state), state.n)

function molar_volumes(state::ThermoState{L}) where L 
    x = state.n[:]./sum(state.n[:])

    mixedvol = volume(state.model.mixed, state.P, state.T, x, phase=state.phase)
    purevol  = volume.(state.model.pure, state.P, state.T, 1.0, phase=state.phase)
    return Species{L}(purevol.*(mixedvol/sum(purevol.*x)))
end
molar_volumes(::Type{<:Species{L}}, state::ThermoState) where L = molaravgs(Species{L}, molar_volumes(state), state.n)

function readvalues(d::Dict{<:Any,<:ET}, obj::ThermoState{L}) where {L,ET}
    getter = Base.Fix1(getindex,d)
    return ThermoState{L,ET}(
        model = obj.model,
        T = getter(obj.T),
        P = getter(obj.P),
        n = Species{L,ET}(getter.(obj.n[:])),
        phase = obj.phase
    )
end



#=======================================================================================
# Test code
=======================================================================================#
#=
clapmap = Dict{Symbol,String}(
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

L = Tuple(collect(keys(clapmap)))
model = ThermoModel{L}(clapmap)
state = ThermoState{L,Float64}(model=model, T=273.15, P=101.3, n=Species{L}(rand(length(L))), phase=:gas)

molar_weights(model)
molar_volumes(state)
=#

