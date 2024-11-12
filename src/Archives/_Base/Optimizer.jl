using LinearAlgebra

Base.@kwdef struct ValGradHess{T<:Number}
    value    :: T
    gradient :: Vector{T}
    hessian  :: Symmetric{Float64,Array{Float64,2}}
end


function box_constrained_min(f_y∇H::Function, x0; f_converged = (Δx, Δy)-> false , maxSteps=10, LB=nothing, UB=nothing)
    elType = eltype(x0)
    n = length(x0)
    y∇H = f_y∇H(x0) #Second order derivatives
    y0  = y∇H.value

    stepN = 0
    converged = false

    LBx = something(LB, fill(-Inf,n))
    UBx = something(UB, fill(Inf,n))

    while (stepN <= maxSteps) && (!converged)
        Δx = box_newt_Δx(x0, y∇H, LBx, UBx)
        x = x0 + Δx

        y∇H = f_y∇H(x)
        Δy = y∇H.value - y0

        converged = f_converged(Δx, Δy)
        y0=y∇H.value
        x0=x
        stepN += 1
    end

    return (x=x0, y∇H=y∇H)
end

function box_newt_Δx(x0, y∇H, LB, UB)
    (∇,H) = (y∇H.gradient, y∇H.hessian)

    #Take pure Newton step (Δx = (x-x0)) by solving for Δx
    # y  ≈ y0 + Δx*∇ + 0.5*Δx'*H*Δx;
    # ∇y ≈ ∇ + H*Δx
    # Δx = -H\∇

    ΔLB = LB.-x0
    ΔUB = UB.-x0

    #return any indices where gradient isn't zero and doesn't collide with the constraint
    function isfree(Δx)
        s = -(∇ .+ H*Δx) #gradient step direction
        isFree = trues(length(s))
        isFree .= isFree .&  .!( (s .< 0) .& isapprox.(Δx, ΔLB) ) #Check lower bound
        isFree .= isFree .&  .!( (s .> 0) .& isapprox.(Δx, ΔUB) ) #Check upper bound
        return isFree
    end
   
    Δx  = clamp.( -(H \ ∇), ΔLB, ΔUB )
    a = isfree(Δx)

    #Determine which variables are free (a) = not yet optimized
    nIter = 0
    searching = any(a)

    while searching && (nIter < 100) #Clamp the inactive states
        c = .! a
        Δxa = -( H[a,a] \ (∇[a] + H[a,c]*Δx[c]) )
        Δx[a] .= clamp.(Δxa, ΔLB[a], ΔUB[a])

        nIter += 1
        aNew = isfree(Δx)
        searching = !(a == aNew)
        a = aNew
    end
    println("$(nIter) active set iterations")

    return Δx
end


#=
using ForwardDiff
using Statistics
n = 3000
p = 500
J = randn(p,n)
R = Diagonal(ones(p))
P = Diagonal(ones(n))
Si = inv( J'*(R\J) + P)

f_objective(x) = x'*Si*x
LB = fill(-Inf, n)
UB = fill(Inf, n)

LB .= 0.1*randn(n)
UB .= 10
function f_y∇H(x)
    return ValGradHess{eltype(x)}(
        value = f_objective(x),
        gradient = 2*Si*x,
        hessian = Symmetric(2*Si)
    )
end

x0 = LB .+ 0.1
@time xn = box_constrained_min(f_y∇H, x0, LB=LB, UB=UB, maxSteps=3)



using Optim
function g_objective!(G,x)
    G .= 2*Si*x
    return G
end

@time xn2 = optimize(f_objective, g_objective!, LB, UB, x0, Fminbox(LBFGS()))
=#