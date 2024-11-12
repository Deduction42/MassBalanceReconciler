
using ForwardDiff
using FiniteDiff
using Statistics

f_nL(K, n, RL) = n./(K./RL .+ 1)

function f_RL(K, n, RL)
    nL = f_nL(K, n, RL)
    return 1/(sum(n)/sum(nL) - 1)
end

function f_RL(K, n)
    RL = eltype(K+n)(1) #High initial guess (function is nearly linear)
    err(x) = x - f_RL(K, n, x)

    for ii in 1:10 #n iterations of Newton's method, it's nearly a straight line at high values
        Δ  = err(RL)/ForwardDiff.derivative(err, RL)
        if !isfinite(Δ)
            break
        end
        RL = abs(RL - Δ)
    end

    return RL
end

n = [10, 20, 30, 10]
K = [10, 0.01, 2.1, 1.5]

RL = f_RL(K,n)

f_nL(K,n) = f_nL(K,n, f_RL(K,n))

@time g1 = FiniteDiff.finite_difference_gradient(x->f_RL(x,n), K)
@time g2 = ForwardDiff.gradient(x->f_RL(x,n), K)

nL = f_nL(K,n, f_RL(K, n))
nG = n.- nL

nG./nL * (sum(nL)/sum(nG))

#The methodology for K values can be found here
# http://www.jmcampbell.com/tip-of-the-month/2006/09/how-to-determine-k-values/
# Basically assume Temperature and Pressure and get K-Values

using PyPlot; pygui(true)
x = LinRange(-2,5,300)
plot(x,err.(x))