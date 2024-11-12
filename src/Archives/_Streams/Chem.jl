using DataFrames
using DelimitedFiles

# ====================================================================================
# General polynomial correlations (frequently used for thermodynamics)
# ====================================================================================
Base.@kwdef struct Polynomial{N}
    Powers :: NTuple{N,Int32}
    Params :: NTuple{N,Float64}
    Domain :: NTuple{2,Float32}
end

function poly_fx(x0::Number, p::Polynomial)
    y = zero(typeof(x))
    x = clamp(x0, p.Domain...)
    for (a, k) in zip(p.Params, p.Powers)
        if !iszero(a)
            y += a*x^k
        end
    end

    return y
end

function poly_∫fdx(x0::NTuple{2, <:Number}, p::Polynomial)
    y = zero(typeof(x[1]))
    x = Tuple(clamp(xi, p.Domain...) for xi in x0) #No extrapolation

    for (a, k) in zip(p.Params, p.Powers)
        if !iszero(a)
            if k == (-1) #Special integral case of integrating 1/x
                y += a*( log(x[2]) - log(x[1]) )
            else
                ki = k+1 #Integral exponent
                y += a/ki*(x[2]^ki - x[1]^ki)
            end
        end
    end

    #Constant integration for extrapolation
    if x0[1] != x[1]
        y += poly_fx(x[1],p)*(x[1]-x0[1])
    end

    if x0[2] != x[2]
        y += poly_fx(x[2],p)*(x0[2]-x[2])
    end

    return y
end



# ====================================================================================
# Chemical Data
# ====================================================================================

#defines a "get_component_indexer" function that returns a constant; changes cause recompilation
function build_component_indexer(names::NTuple{N, Symbol}) where N
    num = Tuple(1:N)
    ex = :(get_component_indexer() = NamedTuple{Tuple($names)}($num))
    eval(ex)
end

Base.@kwdef mutable struct ChemicalData
    Names :: Vector{Symbol}  = []
    R     :: Float64 = 8.31446261815324
    MolWeight  :: Vector{Float64} = []
    LiqVol :: Vector{Float64} = []
    SolVol :: Vector{Float64} = []
    MeltTemp :: Vector{Float64} = []
    CritTemp :: Vector{Float64} = []
    ΔH_fus   :: Vector{Float64} = []
    ΔH_vap :: Vector{Polynomial{4}} = []
    Cp_Gas :: Vector{Polynomial{5}} = []
    Cp_Liq :: Vector{Polynomial{5}} = []
    Cp_Sol :: Vector{Polynomial{1}} = []
end

function compile_chemical_data()
    thisdir(str::String) = joinpath(@__DIR__ , str)
    dfBasic = convert_columns!(quick_read_csv(thisdir("Chem_Basic.csv")))

    polyType = eltype(fieldtype(ChemicalData, :Cp_Gas))
    Cp_Gas = build_polynomial(thisdir("Chem_CpGas.csv"), dfBasic.Name, polyType)

    polyType = eltype(fieldtype(ChemicalData, :Cp_Liq))
    Cp_Liq = build_polynomial(thisdir("Chem_CpLiq.csv"), dfBasic.Name, polyType)

    polyType = eltype(fieldtype(ChemicalData, :Cp_Sol))
    Cp_Sol = build_polynomial(thisdir("Chem_CpSol.csv"), dfBasic.Name, polyType)

    return ChemicalData(
        Names       = dfBasic.Name,
        MolWeight   = dfBasic[!,"Molecular Weight"],
        LiqVol      = dfBasic[!,"Liquid Vol"],
        SolVol      = dfBasic[!,"Solid Vol"],
        MeltTemp    = dfBasic[!,"Melting Point"],
        CritTemp    = dfBasic[!,"Critical Temperature"],
        Cp_Gas      = Cp_Gas,
        Cp_Liq      = Cp_Liq,
        Cp_Sol      = Cp_Sol
    )
end



# ====================================================================================
# Helper functions
# ====================================================================================

#fileName = joinpath(@__DIR__, "Chem_CpGas.csv")
function build_polynomial(fileName::String, componentNames::AbstractVector, polyType::DataType=Polynomial{5})
    function row2domaintuple(T::DataType, r::DataFrameRow)
        return (T(r[2]), T(r[3]))
    end

    function row2paramtuple(T::DataType, r::DataFrameRow) 
        v = zeros(T, polyType.parameters[1])
        for (ii, x) in enumerate(r[4:end])
            v[ii] = x
        end
        return Tuple(v)
    end

    df = convert_columns!( quick_read_csv(fileName) )
    powers = row2paramtuple(Int32, df[1,:])
    scales = row2paramtuple(Float64, df[2,:])
    params = order_by_column!(df[3:end,:], :Name, componentNames)
    
    N = length(powers)
    polyParams = Polynomial[]
    for r in eachrow(params)
        par = polyType(Powers=powers, Params=row2paramtuple(Float32,r).*scales, Domain=row2domaintuple(Float32,r))
        push!(polyParams, par)
    end

    return polyParams
end



function quick_read_csv(fileName::String)
    rawData = readdlm(fileName, ',', String)
    return DataFrame(rawData[2:end,:], Symbol.(rawData[1,:]))
end

function convert_columns!(df, conversions=Pair[]; defaultConversion=any2float)
    colDict = Dict(conversions)
    colDict["Name"] = x->Symbol(x)

    for colName in names(df)
        fconvert = get(colDict, colName, defaultConversion)
        df[!,colName] = fconvert.(df[!,colName])
    end
    return df
end

function order_by_column!(df, colName, colTarget::Vector)
    indexer  = build_indexer(colTarget)
    orderVec = [indexer[n] for n in df[!,colName]]

    for c in eachcol(df)
        c[orderVec] .= c
    end

    return df
end

function build_indexer(targetVec)
    pairGen = ( p[2]=>p[1] for p in enumerate(targetVec) )
    return (; pairGen...)
end

massvec2molevec(vm::AbstractVector{<:Real}) = vm ./ MOD_ENV.ComponentMW
molevec2massvec(vn::AbstractVector{<:Real}) = vn .* MOD_ENV.ComponentMW

function massfrac2molefrac(vm::AbstractVector{<:Real})
    moles = massvec2molevec(vm)
    return moles / sum(filter(!isnan, moles))
end

function molefrac2massfrac(vn::AbstractVector{<:Real})
    masses = molevec2massvec(vn)
    return masses / sum(filter(!isnan, masses))
end


#C = [0.29,0.08619,1.7016,0.00103,909.8]
#=
C = [1.3554, 4.431, 1.6356, 3.054, 746.4]
C = C.*[1e-5,1e-5,1e-3,1e-5,1]  

f(T) = C[1] + C[2]*( (C[3]/T)/sinh(C[3]/T) )^2 + C[4]*( (C[5]/T) /cosh(C[5]/T) )^2 
vT = 0:10:2000

using PyPlot; pygui(true)
using QuadGK

N = 100000

@time for ii in 1:N
    y = f(369.0)
end

@time for ii in 1:N
    (x,w) = gauss(6, 200, 500)    
    y = sum(f.(x) .* w)
end


plot(vT, [f.(vT) X*β])
=#
