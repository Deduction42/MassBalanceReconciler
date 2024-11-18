include("_ThermoModel.jl")

const GhgSpecies{T} = Species{(:CO2, :CH4, :N2O, :other), T, 4} where T

const GHG_MAP = GhgSpecies{Union{Symbol,Nothing}}(
    CO2   = :co2,
    CH4   = :methane,
    N2O   = nothing,
    other = nothing
)

function GhgSpecies(mixture::Species{L,T}, ghgmap::GhgSpecies{Union{Symbol,Nothing}}) where {L,T}
    function get_ghg(ghg::Symbol)
        translation = ghgmap[ghg]
        if isnothing(translation)
            return zero(T)
        else
            return mixture[translation]
        end
    end

    gases = species(GhgSpecies)[begin:(end-1)]
    gasvals = get_ghg.(gases)

    return GhgSpecies{T}((gasvals..., sum(mixture)-sum(gasvals)))
end

#=======================================================================================
# Test code
=======================================================================================#
model = ThermoModel(collect(keys(CLAPEYRON_MAP)))
L = species(model)
z = Species{L}(rand(length(L)))
ghg = GhgSpecies(z, GHG_MAP)
