
struct Options{T<:AbstractFloat, U<:AbstractVector}

    # Display Options
    display::Bool 
    displayInterval::UInt16

    # Function tolerance
    funcTol::T

    # Check function val for NaN and Inf 
    funValCheck::Bool

    # Initial bounds
    iUB::U
    iLB::U

    # Max Iterations
    maxIters::UInt32

    # Max Stall Iterations 
    maxStallIters::UInt32

    # Max Stall Time 
    maxStallTime::T

    # Max Time
    maxTime::T

    # Objective Limit 
    objLimit::T

    # Use parallel
    useParallel::Bool 

    # Is vectorized 
    #isVectorized::Bool # Not implemented
end

function Options(;display=true, displayInterval=1, funcTol=1e-6,
    funValCheck=true, iUB::V=nothing, iLB::V=nothing, maxIters=1000, maxStallIters=25,
    maxStallTime=500, maxTime=1800, objLimit=-Inf, useParallel=false) where 
        {V<:Union{Nothing,AbstractArray}}

    T = typeof(funcTol)
    if iUB === nothing
        iUB = Vector{T}([])
    end
    if iLB === nothing
        iLB = Vector{T}([])
    end

    return Options{T,typeof(iUB)}(display, displayInterval, funcTol,
        funValCheck, iUB, iLB, maxIters, maxStallIters, maxStallTime, 
        maxTime, objLimit, useParallel)
end
