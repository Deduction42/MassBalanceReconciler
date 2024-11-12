
function include_folder(folderName::String)
    folderPath = joinpath(@__DIR__, folderName)
    loadFile = joinpath(folderPath, "__load.jl")

    if isfile(loadFile)
        include(loadFile)
    else
        juliaFiles = filter(x->endswith(x,".jl"), readdir(folderPath))
        for fileName in sort(juliaFiles)
            include( joinpath(folderPath,fileName) )
        end
    end
    return nothing
end

include_folder("_Base")
include_folder("_Streams")
include_folder("_Measurements")
include_folder("HeatExchanger")

Base.@kwdef mutable struct EnvironmentType
    ΔT :: Float64 = 30
end

const MOD_ENV = EnvironmentType()


#Generated functions are used to improve perfomance 
#However they need to have all used functions declared beforehand (otherwise this results in a "new world" problem)
#Include generated indexers after all key_states functions have been built
include(joinpath(@__DIR__, "_Base", "GeneratedIndexers.jl"))

