struct Problem{fType,S}

    f::fType  # Objective function
    
    LB::Vector{S}
    UB::Vector{S}

end

# If LB and UB are not specified, use default of +/- 1000
function Problem(objFunc, numVars)

    # Check that nDims > 0
    numVars > 0 || throw(ArgumentError("N must be greater than zero."))

    # Get types
    fType = typeof(objFunc)

    # Initialize lower and upper bounds
    LB = Vector{Int}(undef, numVars)
    UB = Vector{Int}(undef, numVars)
    @inbounds for i in 1:numVars
        LB[i] = -1000
        UB[i] = 1000
    end

    return Problem{fType,Int}(objFunc, LB, UB)
end

function Problem{S}(objFunc, LB, UB) where {S}

    # Check that LB and UB are the same length
    N = length(LB)
    if length(UB) != N
        throw(ArgumentError("Lower bounds and upper bound vectors must be the same length."))
    end

    # Check that LB[i] < UB[i] ∀ i
    @inbounds for i in 1:N
        LB[i] < UB[i] || throw(ArgumentError("LB[i] must be greater than UB[i] for all i."))
    end

    # Get f type
    fType = typeof(objFunc)

    # Initialize lower and upper bounds
    sLB = Vector{S}(LB)
    sUB = Vector{S}(UB)

    return Problem{fType,S}(objFunc, sLB, sUB)
end

function Problem(objFunc, LB::AbstractArray{S}, UB::AbstractArray{S}) where {S}
    return Problem{S,N}(objFunc, LB, UB)
end
