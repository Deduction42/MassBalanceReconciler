
firstFiles = ["Base.jl", "CanonicalTypes.jl"]
lastFiles = setdiff( readdir(@__DIR__), ["__load.jl"; firstFiles] )

for fileName in vcat(firstFiles, lastFiles) 
    include(joinpath.(@__DIR__, fileName))
end