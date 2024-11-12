#https://en.wikipedia.org/wiki/Flash_evaporation
using Roots
using ForwardDiff
using FiniteDiff

function evaporated_frac(n::AbstractArray{<:Number}, K::AbstractArray{<:Number})
    RetType = promote_type( eltype(n), eltype(K) )
    (kMin, kMax) = extrema(K)

    if kMin > 1 #Everything is a gas
        return RetType(1.0)
    elseif kMax < 1 #Everything is a liquid
        return RetType(0.0)
    end

    #Acquire mole fractions
    z  = n./sum(n)
    dK = K .- 1.0

    #Set up Rachford-Rice equation
    function err(β)
        errGen = ( zi*dKi/(1 + β*dKi) for (zi, dKi) in zip(z,dK) )
        return sum( errGen )
    end

    #Set up initial guess and limits for Brent's method
    ϵ = 1e-6
    βLim = sort( 1 ./ (1 .- [kMax, kMin]) )

    β = find_zero(err, βLim+[ϵ,-ϵ], Roots.Brent())
    return clamp(β, RetType(0.0), RetType(1.0))
end

n = 50 .* rand(10)
K = 2 .* rand(10)

f_K(K) = evaporated_frac(n,K)

@time β  = evaporated_frac(n, K)
@time g1 = FiniteDiff.finite_difference_gradient(f_K, K)
@time g2 = ForwardDiff.gradient(f_K, K)

#using PyPlot; pygui(true)
#x = LinRange(0,1,1000)
#plot(x, err.(x))