# !!! This solves the system base don RG, but RL might be a better method
# RG can never be 0, but RL may be (at zero pressure everything is gas)
# Try writing system to be solved in terms of RL and see if we can set the initial guess as RL=0


using ForwardDiff
using FiniteDiff
using Statistics

f_nG(K, n, RG) = n./(1.0./(RG*K) .+ 1)

function f_RG(K, n, RG)
    nG = f_nG(K, n, RG)
    return 1/( sum(n)/sum(nG) - 1)
end

function f_RG(K, n)
    RG = eltype(K+n)(1) #Equal distribution as initial guess
    err(x) = x - f_RG(K, n, x)

    for ii in 1:10 #n iterations of Newton's method, it's nearly a straight line at high values
        Δ  = err(RG)/ForwardDiff.derivative(err, RG)
        if !isfinite(Δ)
            break
        end
        RG = abs(RG - Δ)
    end

    return RG
end

n = [1, 20, 30, 10]
K = [0.1, 0.01, 0.1, 0.1]

RG = f_RG(K,n)

f_nG(K,n) = f_nG(K,n, f_RG(K,n))

@time g1 = FiniteDiff.finite_difference_gradient(x->f_RG(x,n), K)
@time g2 = ForwardDiff.gradient(x->f_RG(x,n), K)

nG = f_nG(K,n, f_RG(K, n))
nL = n.- nG

nG./nL * (sum(nL)/sum(nG))

#The methodology for K values can be found here
# http://www.jmcampbell.com/tip-of-the-month/2006/09/how-to-determine-k-values/
# Basically assume Temperature and Pressure and get K-Values

using PyPlot; pygui(true)
x = LinRange(-2,5,300)
plot(x,err.(x))

