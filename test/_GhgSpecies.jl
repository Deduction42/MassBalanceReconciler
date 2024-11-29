GHG = (:CO2, :CH4, :N2O, :other)

const GhgSpecies{T} = Species{GHG, T, 4} where T

"""
buildmap(mapping::GhgSpecies{Symbol}, species::NTuple{N,Symbol}) where N

Uses "mapping" (which maps GHG to symbols in "species"). 
    Returns a GhgSpecies{Vecotr{Int}} that can be used to quickly look up components of "species"
"""
function buildmap(mapping::GhgSpecies{Symbol}, species::NTuple{N,Symbol}) where N
    findspecies(gas::Symbol) = findall(x->(x==mapping[gas]), species)

    gasmap = GhgSpecies{Vector{Int}}(
        CO2 = findspecies(:CO2),
        CH4 = findspecies(:CH4),
        N2O = findspecies(:N2O),
        other = Vector(1:N)
    )
    setdiff!(gasmap[:other], union(gasmap[:CO2], gasmap[:CH4], gasmap[:N2O]))

    return gasmap
end

"""
totals(mapping::GhgSpecies{<:AbstractVector{<:Integer}}, mixture::Species{L,T}) where {L,T}

Aggegates "species" into GHG categories as a total
"""
function ghgtotals(ghgmap::GhgSpecies{<:AbstractVector{<:Integer}}, mixture::Species{L,T}) where {L,T}
    function get_total(inds)
        return sum(ind->mixture[ind], inds, init=zero(T))
    end
    return GhgSpecies{T}(get_total.(ghgmap))
end

"""
molaravgs(mapping::GhgSpecies{<:AbstractVector{<:Integer}}, mixture::Species{L,T}) where {L,T}

Aggegates specific properies of "species" as GHG categories (using mole fractions)
"""
function ghgmolaravgs(ghgmap::GhgSpecies{<:AbstractVector{<:Integer}}, mixture::Species{L,T}, moles::Species{L}) where {L,T}
    RT = promote_type(T, Float64)
    function get_molaravgs(inds)
        moltot = sum(ind->moles[ind], inds)
        if iszero(moltot)
            return sum(ind->mixture[ind], inds)/length(inds)
        else
            return sum(ind->mixture[ind]*moles[ind]/moltot, inds)
        end
    end
    return GhgSpecies{RT}(get_molaravgs.(ghgmap))
end

#=======================================================================================
# Test code
=======================================================================================#
#=

const GHG_MAP = GhgSpecies{Union{Symbol,Nothing}}(
    CO2   = :co2,
    CH4   = :methane,
    N2O   = nothing,
    other = nothing
)

model = ThermoModel(collect(keys(CLAPEYRON_MAP)))
L = species(model)
z = Species{L}(rand(length(L)))
ghg = GhgSpecies(z, GHG_MAP)
=#