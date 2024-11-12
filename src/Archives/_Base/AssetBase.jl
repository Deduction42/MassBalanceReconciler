# ==============================================================================
# General functions
# ==============================================================================
abstract type AbstractAsset{T} end
abstract type AbstractSpecs end
abstract type AbstractProcessAsset{T} <: AbstractAsset{T} end
abstract type AbstractConnector{T} <: AbstractAsset{T} end

struct Bounds
    LB :: Float64
    UB :: Float64
end

Bounds(X::AbstractArray) = Bounds(extrema(X)...)
Bounds(X::Bounds) = X

mutable struct GenericAsset{T} <: AbstractAsset{T}
    id :: Symbol
end

## !!! READ THIS!!! It tells you how to write "eltype" correctly
#https://docs.julialang.org/en/v1/manual/methods/


import Base.eltype
basetype(T::Type) = typename(T).wrapper
basetype(x::AbstractAsset) = basetype(typeof(x))

typename(T::DataType) = T.name
typename(T::UnionAll) = typename(T.body)

paramtypes(x) = typeof(x).parameters
paramtypes(T::DataType) = T.parameters
paramtypes(T::UnionAll) = example_asset_type(T).parameters
eltype(::Type{<:AbstractAsset{T}}) where {T} = T
eltype(x::AbstractAsset{T}) where {T} = T

any2float(x::String)  = something(tryparse(Float64,x), NaN64)
any2float(x::Number)  = Float64(x)
any2float(x::Any)     = NaN64
any2float(x::AbstractArray) = any2float.(x)

nonzero(x::Number, ϵ=1e-9) = iszero(x) ? ϵ : x
killnan(x::T, default=0) where T<:Number = isnan(x) ? T(default) : x

#Find all roots of an abstract type
function root_types(T::Type)
    vT = subtypes(T)
    if isempty(vT)
        return T
    else
        return vcat( root_types.(vT)... )
    end
end

# ==============================================================================
# Generic prediction ("astype" is set to the actual type as default)
# ==============================================================================
function predict_state!(asset1::ObjType, asset0::ObjType) where ObjType
    predict_state!(asset1, asset0, ObjType)
end

# ==============================================================================
# Inheritance
# ==============================================================================
macro inherit_asset_fields(assetExp::Expr)
    BaseAsset   = Core.eval(@__MODULE__, assetExp.args[1])
    typeSymbols = assetExp.args[2:end]

    typevar(x::Symbol) = TypeVar(x, Union{}, Any)
    AssetType   = BaseAsset{ typevar.(typeSymbols)... }

    return esc(build_fieldlist_expr(AssetType, excluding=[:specs]))
end

macro inherit_asset_fields(assetExp::Symbol)
    AssetType = Core.eval(@__MODULE__, assetExp)
    return esc(build_fieldlist_expr(AssetType, excluding=[:specs]))
end

macro inherit_spec_fields(assetExp::Symbol)
    AssetType = Core.eval(@__MODULE__, assetExp)
    return esc(build_fieldlist_expr(AssetType))
end

function get_type_expression(X::DataType)
    params = Symbol.(X.parameters)
    if isempty(params)
        return Symbol(X.name.wrapper)
    else
        return Expr(:curly, Symbol(X.name.wrapper), params...)
    end
end

function get_type_expression(X::TypeVar)
    return X.name
end

function get_type_expression(X::Symbol)
    return X
end

function build_fieldlist_expr(AssetType::DataType; excluding=Symbol[])
    ex = Expr(:block)
    for fName in filter(x-> !(x in excluding), fieldnames(AssetType) )
        fTypeExpr = get_type_expression(fieldtype(AssetType, fName))
        push!(ex.args, :( $fName :: $fTypeExpr ) )
    end
    return ex
end




# ==============================================================================
# Getting and setting indices to emulate vector behaviour better
# ==============================================================================
import Base.getindex
import Base.setindex!
const IterableType = Union{Array, Tuple, Base.Generator}

valuegen(x) = (x[fn] for fn in propertynames(x))

getindex(Obj::AbstractAsset, fn::Symbol) = getproperty(Obj, fn)
getindex(Obj::AbstractAsset, vfn::IterableType) = (Obj[fn] for fn in vfn)
getindex(Obj::AbstractAsset, vfn::Colon) = getindex(Obj, propertynames(Obj))

setindex!(Obj::AbstractAsset, x, fn::Symbol) = setproperty!(Obj, fn, x)
function setindex!(Obj::AbstractAsset, X::IterableType, vfn::IterableType)
    if length(X) != length(vfn)
        error("Target value indices and assigned values must be of the same length")
    else
        for (x, fn) in zip(X, vfn)
            Obj[fn] = x
        end
    end
    return nothing
end
setindex!(Obj::AbstractAsset, x, vfn::Colon) = setindex!(Obj, x, propertynames(Obj))



# ==============================================================================
# Build a similar asset with a new default value for each state
# ==============================================================================
function fill_asset(x, oldAsset::AbstractAsset)
    #n  = length(MOD_ENV.ComponentNames)
    passThroughs = passthrough_fields(typeof(oldAsset))
    oldEltype = eltype(oldAsset)
    BaseType = basetype(oldAsset)
    
    function fill_field(fn, fval)
        if fn in passThroughs
            return fval
        elseif fval isa oldEltype
            return x
        elseif fval isa Array{oldEltype}
            return fill(x, size(fval))
        elseif fval isa AbstractAsset
            return fill_asset(x, fval)
        else
            return fval
        end
    end

    gen = ( fill_field(fn, oldAsset[fn]) for fn in fieldnames(BaseType) )
    return BaseType(gen...)
end

# ==============================================================================
# Build out a generic DataType example for a (possibly UnionAll) datatype
# ==============================================================================
function example_asset_type(::Type{assetType}, ElType::DataType=Float64) where (assetType <: AbstractAsset)
    #Extract the base type
    BaseType = basetype(assetType)

    #Most assets only have the eltype as a parameter
    TN = BaseType{ElType}
    if TN isa DataType
        return TN
    end

    #Every parameter after the eltype is an asset
    #Keep adding assets to the abstract type (UnionAll) until it is fully parameterized (DataType)
    while true 
        Ti = TN{GenericAsset{ElType}}
        if Ti isa DataType
            return Ti
        else
            Tn = Ti
        end
    end
end

# ==============================================================================
# Determine which components are in the BOM, and which ones are inner
# ==============================================================================
function all_components(::Type{T}) where T <: AbstractAsset
    exampleType = example_asset_type(basetype(T))
    return fields_of_type(exampleType, AbstractAsset)
end

function all_components(X::AbstractAsset)
    return all_components(typeof(X))
end

#Asset components that you DONT want showing up in the BOM (default is nothing)
function inner_components(X::AbstractAsset)
    return inner_components(typeof(X))
end

#Asset components that you DONT want showing up in the BOM (default is nothing)
function inner_components(::Type{T}) where T<: AbstractAsset
    return Symbol[]
end


#Asset components that you DO want to show up in the BOM (usually streams)
function bom_components(::Type{T}) where T<: AbstractAsset
    allComponents   = all_components(T)
    innerComponents = inner_components(T)

    if isempty(innerComponents)
        return allComponents
    else
        gen = ( fn for fn in allComponents if !(fn in innerComponents) ) 
        return Tuple(gen)
    end
end

function bom_components(X :: AbstractAsset)
    return bom_components(typeof(X))
end
