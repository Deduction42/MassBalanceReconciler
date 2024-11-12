using Optim
using ForwardDiff

mutable struct BomEnvType
    Archives :: Dict{DataType, Any}
    Indexer  :: NamedTuple
end



# ==============================================================================
# Gives us the ability to reuse a BOM instead of allocating a new one
# BOM Archives is a dict containing various BOM types
# BOM Must have an Int32 entry
# ==============================================================================
function reuse_bom_type!(BOM_Archives, ::Type{T}) where T <: Type
    if !haskey(BOM_Archives, T)
        indexBom = BOM_Archives[Int32]
        BOM = (t0=new_eltype(T, indexBom), t1=new_eltype(T, indexBom))
        BOM_Type_Archive[T] = BOM
        return BOM
    else
        return BOM_Achives[T]
    end
end

# ==============================================================================
# BOM prediction set
# ==============================================================================
function predict_bom_state(stateVec::AbstractVector, indexer::NamedTuple, BOM_Archives::Dict)
    bom = reuse_bom_type( BOM_Archives, eltype(stateVec) )
    set_bom_state!(bom.t0, stateVec, indexer)

    #Filter out connectors from the equipment index list
    equipmentInd = [ ii for ii in eachindex(bom.t0) if !(bom.t0[ii] isa AbstractConnector) ]

    Threads.@threads for ii in equipmentInd
        predict_state!(bom.t1[ii], bom.t0[ii])
    end

    return bom.t1
end

# ==============================================================================
# Updating Assets and the BOM
# ==============================================================================
function set_bom_state!(BOM::NamedTuple, stateVec::AbstractVector, indexer::NamedTuple)
    for ii in 1:length(indexer)
        set_state!(BOM[ii], stateVec, indexer[ii])
    end

    return BOM
end

function get_bom_state!(stateVec::AbstractVector, BOM::NamedTuple, indexer::NamedTuple)
    for name in keys(BOM)
        get_state!(stateVec, BOM[name], indexer[name])
    end

    return stateVec
end

# Make a BOM of indexers from a BOM in either dictionary or NamedTuple format
function bom_indexer(BOM::Union{<:AbstractDict,NamedTuple})
    lastInd = Index_Type(0)
    bomKeys = keys(BOM)
    bomVals = (build_indexer!(lastInd, Asset) for Asset in values(BOM))

    bomIND  = NamedTuple{Tuple(bomKeys)}(Tuple(bomVals))
    return (bom = connect_bom!(bomIND), n = lastInd.Val) 
end

#Change element types in a BOM
function new_eltype(elType::Type, oldBom::NamedTuple)
    newVals = ( new_eltype(elType, Asset) for Asset in values(oldBom) )
    newBom  = NamedTuple{propertynames(oldBom)}( Tuple(newVals) )
    connect_bom!(newBom)
    return newBom
end

# ==============================================================================
# Connect the bom so that the components point to the bom indices
# ==============================================================================
function connect_bom!(bom::NamedTuple)
    for asset in bom
        connect_asset!(asset, bom)
    end
    return bom
end

function connect_asset!(asset::AbstractAsset, bom::NamedTuple)
    for fn in bom_components(typeof(asset))
        asset[fn] = bom[asset[fn].id]
    end
    return bom
end






if false

    S1_Str = ProcessStream{String}()
    S1_Str.id = Symbol("S1")
    S1_Str.type = ProcessStream_Water

    LastInd = Index_Type(0)
    S1_Ind  = new_eltype(Int32, S1_Str)
    S1_Ind  = fill_indexer!(LastInd, S1_Ind)

    Indexer = S1_Ind
    StateVec = [1.0,2,3,4,5,6]

    BOM_Str = Dict(:S1=>S1_Str)
    BOM_Ind = bom_indexer(BOM_Str)

    BOM_Float = new_eltype(Float64, BOM_Ind)
    set_bom_state!(BOM_Float, StateVec, BOM_Ind)
    BOM_Float[:S1]

    using Optim

    function f_obj(x)
        BOM = reuse_bom_type(eltype(x), BOM_Ind)
        set_bom_state!(BOM, x, BOM_Ind)

        ε = sum(x.-StateVec)^2
        ε += (BOM[:S1].temperature - 25)^2 / 0.1
        ε += (BOM[:S1].pressure - 150)^2 / 5
        ε += (sum(BOM[:S1].mass_flows) - 100)^2 / 10

        return ε
    end

    function g_obj!(g, x)
        g .= ForwardDiff.gradient(f_obj, x)
    end

    g_obj!(StateVec*0, StateVec*0)
    f_obj(StateVec*0)

    Result = optimize(f_obj, g_obj!, StateVec*0, LBFGS())
    StateVec = Result.minimizer
    set_bom_state!(BOM_Float, StateVec, BOM_Ind)
    BOM_Float[:S1]
end
