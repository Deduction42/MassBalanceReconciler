using Optim
using Zygote

function errorgradient(statevec::AbstractVector, plant::PlantState)
    return gradient(x->negloglik(x, plant), statevec)[1]
end

