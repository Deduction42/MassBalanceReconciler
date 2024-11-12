Base.@kwdef struct AssetComponents
    BaseType :: UnionAll
    BOM   :: Vector{Symbol}
    Inner :: Vector{Symbol}
end

# ==============================================================================
# Builds a BOM from a list of dictionaries
# ==============================================================================
function build_bom(dictList::Vector; elType=Float64)
    keyList  = [Symbol(dictItem["id"]) for dictItem in dictList]
    assetBOM = Dict{Symbol,Any}( Pair.(keyList, dictList) )

    for asset in values(assetBOM)    
        build_asset!(assetBOM, asset, elType=elType)
    end

    #Build NamedTuple BOM with assets sorted in order (for easy inspection)
    pairGen = ( k=>assetBOM[k] for k in sort(keyList) )

    return (; pairGen...)
end


# ==============================================================================
# Build an asset from a dictionary
# ==============================================================================
function build_asset!(assetBOM::Dict{Symbol}, id::Symbol; elType=Float64)
    return build_asset!(assetBOM, assetBOM[id]; elType=elType)
end

function build_asset!(assetBOM::Dict{Symbol}, id::String; elType=Float64)
    return build_asset!(assetBOM, assetBOM[Symbol(id)]; elType=elType)
end

function build_asset!(assetBOM::Dict{Symbol}, assetDict::Dict{String}; elType=Float64)
    assetType = eval(Symbol(assetDict["class"]))
    asset = build_asset(assetType, assetDict, BOM=assetBOM, elType=elType)
    assetBOM[asset.id] = asset
    return asset
end

function build_asset!(assetBOM::Dict{Symbol}, asset::AbstractAsset; elType=Float64)
    return asset
end

function build_asset(::Type{T}, assetDict::Dict{String}; BOM=nothing, elType=Float64) where T <: AbstractAsset
    comp = asset_components(T)
    baseType = comp.BaseType
    exampleType = example_asset_type(baseType, elType)

    n = length(CHEM.Names)

    function get_with_default(assetDict, fs)
        try 
            return assetDict[fs]
        catch
            return default_value(exampleType, Symbol(fs))
        end
    end

    function build_field(fn)
        fs = String(fn)
        fT = fieldtype(exampleType, fn)

        if fT <: AbstractAsset 
            if isnothing(BOM) || (fn in comp.Inner)
                return build_asset(fT, assetDict[fs], BOM=nothing, elType=elType)
            else
                return build_asset!(BOM, assetDict[fs], elType=elType)
            end

        elseif fT <: AbstractSpecs
            return build_specs(fT, get(assetDict, fs, nothing))

        elseif fT <: Vector{<:elType}
            return elType.(get_with_default(assetDict, fs))

        elseif fT <: elType
            return elType(get_with_default(assetDict, fs))

        elseif fT <: DataType
            return eval(Symbol(assetDict[fs]))

        else
            return fT(get_with_default(assetDict, fs))   
        end
    end

    #Build the asset and assign it to the BOM
    assetGen = Tuple( build_field(fn) for fn in fieldnames(T) )
    
    return baseType{elType}( assetGen... )
end

#Builds a specification object from the dictionary
function build_specs(::Type{specT}, objDict::Dict{String}) where specT <: AbstractSpecs
    objGen = ( objDict[String(fn)] for fn in fieldnames(specT) )
    return specT(objGen...)
end

function build_specs(::Type{specT}, objDict::Nothing) where specT <: AbstractSpecs
    error("Couldn't find required specifications for object type \"$(specT)\"")
end

#Get the list of components for an asset
function asset_components(assetDict::Dict)
    return asset_components( assetDict[:class] )
end

function asset_components(typeId::Union{String,Symbol})
    return asset_components( eval(Symbol(typeId)) )
end

function asset_components(AssetType::UnionAll)
    BaseType = basetype(AssetType)
    return AssetComponents(BaseType=BaseType, BOM=bom_components(BaseType), Inner=inner_components(BaseType))
end


function str2sym_dict(D::Dict)
    return Dict( Symbol(p[1])=>str2sym_dict(p[2]) for p in pairs(D) )
end

function str2sym_dict(x::String)
    return Symbol(x)
end

function str2sym_dict(x::Any)
    return x
end
