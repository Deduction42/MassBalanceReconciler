using MassBalanceReconciler
using JSON3
using Revise

plantjson = open(joinpath(@__DIR__,"plant.json")) do fh
    JSON3.read(fh, numbertype=Float64)
end

plantinfo = PlantInfo(plantjson)

plantstate = PlantState(plantinfo)