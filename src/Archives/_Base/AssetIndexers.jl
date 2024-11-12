

# ==============================================================================
# Generic method to fill an indexer
# ==============================================================================
mutable struct Index_Type
    Val :: Int32
end

function iterate_indexer!(LastInd :: Index_Type)
    LastInd.Val += 1
    return LastInd.Val
end


function build_indexer!(lastInd::Index_Type, asset::AbstractAsset)
    indexer = fill_asset(Int32(0), asset)
    fill_indexer!(lastInd, indexer)
    return indexer
end

function fill_indexer!(lastInd::Index_Type, indexer::AbstractAsset)
    #Index inner components
    for fn in inner_components(typeof(indexer)) #Index inner components
        fill_indexer!(lastInd, indexer[fn])
    end

    #Index key states
    fill_indexer_shallow!(lastInd, indexer)

    return indexer
end

function fill_indexer_shallow!(lastInd::Index_Type, indexer::AbstractAsset)
    #Fill in key states (which can be arrays)
    for fn in key_states(typeof(indexer))
        if fieldtype(typeof(indexer), fn) <: AbstractArray
            fArray = getproperty(indexer, fn)
            for ii in eachindex(fArray)
                fArray[ii] = iterate_indexer!(lastInd)
            end

        else
            setproperty!(indexer, fn, iterate_indexer!(lastInd))

        end
    end

    return indexer
end


# ==============================================================================
# Keeping track of indexers
# ==============================================================================
function valid_fieldtypes(aT::UnionAll)
    return NamedTuple{fieldnames(aT)}(fieldtypes(aT))
end

function valid_fieldtypes(aT::DataType)
    bT = basetype(aT)
    return NamedTuple{fieldnames(bT)}(fieldtypes(bT))
end

function valid_fieldtypes(asset::AbstractAsset)
    return valid_fieldtypes(typeof(asset))
end

function fields_of_type(aT::Type, fT::Type)
    return [fn for fn in fieldnames(aT) if fieldtype(aT,fn) <: fT]
end

function fields_of_type(X0::AbstractAsset, fT::Type)
    return fields_of_type(typeof(X0), fT)
end
        
function default_getindex(v::AbstractArray{<:Number}, ind::Integer; default=0)
    if ind == 0
        return eltype(v)(default)
    else
        return v[ind]
    end
end

function default_getindex(v::AbstractArray{<:Number}, ind::AbstractArray{<:Integer}; default=0)
    return [ default_getindex(v, ii, default=default) for ii in ind ]
end

# ==============================================================================
# Fallback methods for key states, compoments and inner components of assets
# ==============================================================================
function key_states(X::AbstractAsset)
    return key_states(typeof(X))
end

function key_states(::Type{assetType}) where assetType <: AbstractAsset
    baseType = basetype(assetType)
    T = Some{:T}
    exampleType = example_asset_type(baseType, T)
    return fields_of_type(exampleType, Union{T, Array{T}})
end



#Fields you automatically pass through when changing eltype
function passthrough_fields(::Type{T}) where T <: AbstractAsset
    ispassthrough(fn) = (fieldtype(T, fn) <: Union{Symbol, AbstractSpecs, Type})
    gen = (fn for fn in fieldnames(T) if ispassthrough(fn))
    return Tuple(gen)
end

