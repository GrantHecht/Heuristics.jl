struct Problem{fType,S,N}

    f::fType  # Objective function
    
    LB::SizedArray{Tuple{N},S,1,1,Vector{S}}
    UB::SizedArray{Tuple{N},S,1,1,Vector{S}}

end

# If LB and UB are not specified, use default of +/- 1000
function Problem{N}(objFunc) where {N}

    # Check that nDims > 0
    N > 0 || throw(ArgumentError("N must be greater than zero."))

    # Get types
    fType = typeof(objFunc)

    # Initialize lower and upper bounds
    LB = SizedArray{Tuple{N},Int,1,1,Vector{Int}}(undef)
    UB = SizedArray{Tuple{N},Int,1,1,Vector{Int}}(undef)
    @inbounds for i in 1:N
        LB[i] = -1000
        UB[i] = 1000
    end

    return Problem{fType,Int,N}(objFunc, LB, UB)
end

function Problem{S,N}(objFunc, LB, UB) where {S,N}

    # Check that N > 0
    N > 0 || throw(ArgumentError("N must be greater than zero."))

    # Check that LB[i] < UB[i] ∀ i
    @inbounds for i in 1:N
        LB[i] < UB[i] || throw(ArgumentError("LB[i] must be greater than UB[i] for all i."))
    end

    # Get f type
    fType = typeof(objFunc)

    # Initialize lower and upper bounds
    sLB = SizedArray{Tuple{N},S,1,1,Vector{S}}(LB)
    sUB = SizedArray{Tuple{N},S,1,1,Vector{S}}(UB)

    return Problem{fType,S,N}(objFunc, sLB, sUB)
end

function Problem{N}(objFunc, LB::AbstractArray{S}, UB::AbstractArray{S}) where {S,N}
    return Problem{S,N}(objFunc, LB, UB)
end
