# ==============================================================================
# Add custom streams
# ==============================================================================
folderName = joinpath(@__DIR__, "Custom Streams")
for fileName in readdir(folderName)
    include(joinpath(folderName,fileName))
end


# ==============================================================================
# Manual dispatch method for streams according to its Substance class
# Generated functions are used in order to dispatch quickly
# Must be rebuilt if any new stream types are added
# ==============================================================================
function _subclass_dispatch_expression(typeList)
    (firstType, remainingTypes) = Iterators.peel(typeList)

    #First type initalizes the if block
    code = :(if x.subclass === $firstType
        callfunc($firstType, x)
    end)

    #The next few iterations append the elseif blocks
    args = code.args
    for nextType in remainingTypes
        clause = :(if x.subclass === $nextType # use `if` so this parses, then change to `elseif`
                    callfunc($nextType, x)
                end)
        clause.head = :elseif
        push!(args, clause)
        args = clause.args #Recursive extention of the this argument's argument
    end

    #The final clause resorts to plain multiple dispatch
    push!(args, :(callfunc(x.subclass,x)))
    return code
end


@generated function stream_dispatch(callfunc, x::ProcessStream)
    typeList = Tuple(root_types(Substance))
    return _subclass_dispatch_expression(typeList)
end

@generated function stream_dispatch(callfunc, x::ProcessStream, ::Type{TypeList}) where TypeList <: Tuple
    typeList = Tuple(TypeList.parameters)
    return _subclass_dispatch_expression(typeList)
end


#=
function stream_dispatch(callfunc, x::ProcessStream)
    return callfunc(x.subclass, x)
end
=#