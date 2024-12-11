using Clapeyron
using LinearAlgebra
using SparseArrays
include("_Species.jl")

import Clapeyron.molecular_weight
import ForwardDiff

abstract type AbstractThermo{L} end
species(::Type{<:AbstractThermo{L}}) where L = L
species(x::AbstractThermo{L}) where L = L

#=======================================================================================
# Thermodynamic models for individual species (including composites)
=======================================================================================#
@kwdef struct ThermoSpecies
    model::PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}
    fracs::Vector{Float64}
end

ThermoSpecies(name::String) = ThermoSpecies(PR(name),[1.0])

function ThermoSpecies(names::AbstractVector{String}, fracs::AbstractVector{<:Real})
    return ThermoSpecies(PR(names), fracs/sum(fracs))
end

molecular_weight(model::ThermoSpecies) = molecular_weight(model.model, model.fracs)

#=======================================================================================
# Thermodynamic models for mixtures
=======================================================================================#

@kwdef struct ThermoModel{L, N} <: AbstractThermo{L}
    pure    :: Species{L, ThermoSpecies, N}
    mixed   :: PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}
    molmap  :: SparseMatrixCSC{Float64, Int64}
    ThermoModel{L,N}(x...) where {L,N} = new{L, length(L)}(x...)
    ThermoModel{L}(x...) where L = new{L, length(L)}(x...)
end

function ThermoModel{L}(substances::AbstractVector{ThermoSpecies}) where L
    mixnames = String[]
    for subst in substances
        union!(mixnames, subst.model.components)
    end
    
    mixed  = PR(mixnames)
    mixmap = Dict{String,Int}(k=>i for (i,k) in enumerate(mixnames))
    molmap = zeros(Float64, (length(mixnames), length(L)))

    for (icol, subst) in enumerate(substances)
        for (comp, frac) in (subst.model.components .=> subst.fracs)
            molmap[mixmap[comp], icol] = frac
        end
    end

    return ThermoModel{L}(Species{L}(substances), mixed, molmap)
end


function ThermoModel{L}(thermomap::Dict{Symbol,String}) where {L}
    species_vec = [thermomap[s] for s in L] 
    models = Species{L}(ThermoSpecies.(species_vec))

    return ThermoModel{L}(models)
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

function Clapeyron.volume(state::ThermoState; T=state.T, P=state.P, n=state.n, phase=state.phase) 
    return volume(state.model.mixed, state.P, state.T, model.molmap*n, phase=phase)
end

function molar_volumes(state::ThermoState{L}) where L 
    moles = state.n[:]
    volfunc(x) = volume(state, n=x)
    vol   = volfunc(moles)
    dvol  = ForwardDiff.gradient(volfunc, moles)
    return Species{L}((vol/dot(dvol, moles)).*dvol)
end

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

claplabels = [
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
]

L = Tuple((p[1] for p in claplabels))
clapmap = Dict(claplabels)

model = ThermoModel{L}(Dict(clapmap))
state = ThermoState{L,Float64}(model=model, T=273.15, P=101.3e3, n=Species{L}(rand(length(L))), phase=:gas)

molar_weights(model)
molar_volumes(state)


