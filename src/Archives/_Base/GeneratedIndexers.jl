# ================================================================================================================
# Indexing generated functions (generated for speed but need to be built last)
# ================================================================================================================
#Sets values of states and components that are not in the BOM (not rigorous, but safe in parallel over the BOM)
@generated function set_state!(state::T1, stateVec::AbstractVector, ind::T2) where {T1 <: AbstractAsset, T2 <: AbstractAsset}
    ex = Expr(:block)

    for fn in key_states(T2)
        if fieldtype(T2, fn) isa Array
            push!(ex.args, :(state.$fn .= stateVec[ind.$fn]) )
        else
            push!(ex.args, :(state.$fn = stateVec[ind.$fn]) )
        end
    end

    for fn in inner_components(T2)
        push!(ex.args, :(set_state!(state.$fn, stateVec, ind.$fn)) )
    end

    return ex
end

#gets the states of the current asset and only inner components (others should be retreived from the BOM)
@generated function get_state!(stateVec::AbstractVector, state::T1, ind::T2) where {T1 <: AbstractAsset, T2 <: AbstractAsset}
    ex = Expr(:block)

    for fn in key_states(T2)
        if fieldtype(T2, fn) isa Array
            push!(ex.args, :(stateVec[ind.$fn] .= state.$fn) )
        else
            push!(ex.args, :(stateVec[ind.$fn] = state.$fn) )
        end
    end

    for fn in inner_components(T2)
        push!(ex.args, :(get_state!(stateVec, state.$fn, ind.$fn)) )
    end

    return ex
end

#Syncrhonizes state1 to take on the same value as state0
@generated function sync_state!(state1::T, state0::T) where {T <: AbstractAsset}
    ex = Expr(:block)

    for fn in key_states(T)
        if fieldtype(T, fn) isa Array
            push!(ex.args, :(state1.$fn .= state0.$fn) )
        else
            push!(ex.args, :(state1.$fn  = state0.$fn) )
        end
    end

    for fn in inner_components(T)
        push!(ex.args, :(sync_state!(state1.$fn, state0.$fn)) )
    end

    return ex
end

#Sets the values of ALL states using the state vector (rigorous but not safe in parallel)
@generated function set_state_deep!(state::T1, stateVec::AbstractVector, ind::T2) where {T1 <: AbstractAsset, T2 <: AbstractAsset}
    ex = Expr(:block)

    for fn in key_states(T2)
        if fieldtype(T2, fn) isa Array
            push!(ex.args, :(state.$fn .= stateVec[ind.$fn]) )
        else
            push!(ex.args, :(state.$fn = stateVec[ind.$fn]) )
        end
    end

    for fn in all_components(T2)
        push!(ex.args, :(set_state_deep!(state.$fn, stateVec, ind.$fn)) )
    end

    return ex
end
