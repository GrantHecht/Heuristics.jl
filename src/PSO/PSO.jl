
struct PSO{T,S,N,M,fType} <: Optimizer

    # Optimization problem
    prob::Problem{fType,S,N}

    # Swarm of particles
    swarm::Swarm{T,N,M}

    # Settings 
    initMethod::Symbol      # Only option now is :Uniform
    updateMethod::Symbol    # Only option now is :MATLAB

    # PSO specific parameters/options
    inertiaRange::Tuple{T,T}
    minNeighborFrac::T 
    selfAdjustWeight::T
    socialAdjustWeight::T 

    # Results data
    results::Results{T,N}

    function PSO{T,N,M}(prob::Problem{fType,S,N}, initMethod::Symbol, 
        updateMethod::Symbol, inertiaRange::Tuple{T,T}, minNeighborFrac::T, 
        selfAdjustWeight::T, socialAdjustWeight::T) where {T,S,N,M,fType}

        # Initialize Swarm 
        swarm = Swarm{T,N,M}(undef)

        # Initialize results
        results = Results{T,N}(undef)

        # Check that init and update methods exist
        if initMethod != :Uniform && initMethod != :LogisticsMap
            throw(ArgumentError("$initMethod is not a PSO initialization method."))
        end
        if updateMethod != :MATLAB 
            throw(ArgumentError("$updateMethod is not a PSO update method."))
        end

        return new{T,S,N,M,fType}(prob, swarm, initMethod, updateMethod, 
                        inertiaRange, minNeighborFrac, selfAdjustWeight, 
                        socialAdjustWeight, results)
    end
end

# ===== Constructors

function PSO{M}(prob::Problem{fType,S,N}; inertiaRange = (0.1, 1.1), minNeighborFrac = 0.25, 
    selfAdjustWeight = 1.49, socialAdjustWeight = 1.49, initMethod::Symbol = :Uniform,
    updateMethod::Symbol = :MATLAB) where {S,N,M,fType}

    # Error checking
    length(inertiaRange) == 2 || throw(ArgumentError("inertiaRange must be of length 2."))
    minNeighborFrac > 0       || throw(ArgumentError("minNeighborFrac must be > 0."))

    # Type info
    T = typeof(1.0)
    nIRange = (T(inertiaRange[1]), T(inertiaRange[2]))

    # Call constructor
    return PSO{T,N,M}(prob, initMethod, updateMethod, nIRange, T(minNeighborFrac),
        T(selfAdjustWeight), T(socialAdjustWeight))
end

# ===== Methods

function _optimize!(pso::PSO, opts::Options)

    # Initialize PSO
    initialize!(pso, opts)

    # Iterate PSO
    iterate!(pso, opts)

    # Display results and return 
    if opts.display
        display(pso.results)
    end

    return pso.results
end

function initialize!(pso::PSO, opts::Options)

    # Initialize position and velocities
    if pso.initMethod == :Uniform
        uniformInitialization!(pso.swarm, pso.prob, opts)
    elseif pso.initMethod == :LogisticsMap
        logisticsMapInitialization!(pso.swarm, pso.prob, opts)
    else
        throw(ArgumentError("PSO initialization method is not implemented."))
    end

    # Evaluate Objective Function
    feval!(pso.swarm, pso.prob.f, opts; init = true)

    # Set swarm global best obj. func. value to Inf 
    pso.swarm.b = Inf

    # Set global best 
    setGlobalBest!(pso.swarm)

    # Initialize neighborhood size
    pso.swarm.n = max(2, floor(length(pso.swarm)*pso.minNeighborFrac))

    # Initial inertia
    if pso.inertiaRange[2] > 0
        pso.swarm.w = (pso.inertiaRange[2] > pso.inertiaRange[1]) ?
                            pso.inertiaRange[2] : pso.inertiaRange[1]
    else
        pso.swarm.w = (pso.inertiaRange[2] < pso.inertiaRange[1]) ?
                            pso.inertiaRange[2] : pso.inertiaRange[1]
    end

    # Initialize self and social adjustment 
    pso.swarm.y₁ = pso.selfAdjustWeight
    pso.swarm.y₂ = pso.socialAdjustWeight

    # Print Status
    if opts.display
        printStatus(pso.swarm, 0.0, 0, 0)
    end

    return nothing
end

function iterate!(pso::PSO, opts::Options)

    # Initialize time and iteration counter
    t0  = time()
    stallT0 = t0
    iters = 0
    stallIters = 0
    fStall = Inf

    # Compute minimum neighborhood size
    minNeighborSize = max(2, floor(length(pso.swarm)*pso.minNeighborFrac))

    # Begin loop
    exitFlag = 0
    while exitFlag == 0

        # Update iteration counter 
        iters += 1

        # Update swarm velocities 
        updateVelocities!(pso.swarm)

        # Update positions 
        step!(pso.swarm)

        # Enforce Bounds
        if !(all(isinf.(pso.prob.LB)) && all(isinf.(pso.prob.UB)))
            enforceBounds!(pso.swarm, pso.prob.LB, pso.prob.UB)
        end

        # Evaluate objective function
        feval!(pso.swarm, pso.prob.f, opts)

        # Update global best
        flag::Bool = setGlobalBest!(pso.swarm)

        # Update stall counter and neighborhood
        if flag
            pso.swarm.c = max(0, pso.swarm.c - 1)
            pso.swarm.n = minNeighborSize
        else 
            pso.swarm.c += 1
            pso.swarm.n = min(pso.swarm.n + minNeighborSize, 
                                length(pso.swarm) - 1)
        end

        # Update inertia
        if pso.swarm.c < 2
            pso.swarm.w *= 2.0
        elseif pso.swarm.c > 5
            pso.swarm.w /= 2.0
        end
        
        # Ensure inertia is in inertia range 
        if  pso.swarm.w > pso.inertiaRange[2]
            pso.swarm.w = pso.inertiaRange[2]
        elseif pso.swarm.w < pso.inertiaRange[1]
            pso.swarm.w = pso.inertiaRange[1]
        end 

        # Track stalling
        if fStall - pso.swarm.b > opts.funcTol
            fStall = pso.swarm.b 
            stallIters = 0
            stallT0 = time()
        else
            stallIters += 1
        end

        # Stopping criteria
        if stallIters >= opts.maxStallIters
            exitFlag = 1
        elseif iters >= opts.maxIters
            exitFlag = 2
        elseif pso.swarm.b <= opts.objLimit
            exitFlag = 3
        elseif time() - stallT0 >= opts.maxStallTime 
            exitFlag = 4
        elseif time() - t0 >= opts.maxTime 
            exitFlag = 5
        end

        # Output Status
        if opts.display && iters % opts.displayInterval == 0
            printStatus(pso.swarm, time() - t0, iters, stallIters)
        end
    end

    # Set results
    setResults!(pso, iters, time() - t0, exitFlag)

    return nothing
end

function setResults!(pso::PSO, iters::Int, Δt::AbstractFloat, exitFlag::Int)
    pso.results.fbest = pso.swarm.b 
    pso.results.xbest .= pso.swarm.d
    pso.results.iters = iters 
    pso.results.time = Δt 
    pso.results.exitFlag = exitFlag

    return nothing
end