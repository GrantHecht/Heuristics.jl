
mutable struct PSO{T,S,fType} <: Optimizer

    # Optimization problem
    prob::Problem{fType,S}

    # Swarm of particles
    swarm::Swarm{T}

    # Settings 
    initMethod::Symbol      # :Uniform or :LogisticsMap
    updateMethod::Symbol    # Only option now is :MATLAB

    # PSO specific parameters/options
    inertiaRange::Tuple{T,T}
    minNeighborFrac::T 
    selfAdjustWeight::T
    socialAdjustWeight::T 

    # Results data
    results::Results{T}

    # Optimizer state flag
    # state = 0 : Not initialized
    # state = 1 : Initialized
    # state = 2 : Optimizing
    # state = 3 : Converged
    state::Int

    # Optimizer time, iteration, and stall parameters
    t0::Float64 
    stallT0::Float64
    iters::Int 
    stallIters::Int
    fStall::T

    function PSO{T}(prob::Problem{fType,S}, numParticles::Integer, initMethod::Symbol, 
        updateMethod::Symbol, inertiaRange::Tuple{T,T}, minNeighborFrac::T, 
        selfAdjustWeight::T, socialAdjustWeight::T) where {T,S,fType}

        # Initialize Swarm 
        N = length(prob.LB)
        swarm = Swarm{T}(N, numParticles)

        # Initialize results
        results = Results{T}(undef, N)

        # Check that init and update methods exist
        if initMethod != :Uniform && initMethod != :LogisticsMap
            throw(ArgumentError("$initMethod is not a PSO initialization method."))
        end
        if updateMethod != :MATLAB 
            throw(ArgumentError("$updateMethod is not a PSO update method."))
        end

        return new{T,S,fType}(prob, swarm, initMethod, updateMethod, 
                        inertiaRange, minNeighborFrac, selfAdjustWeight, 
			socialAdjustWeight, results, 0, 0.0, 0.0, 0, 0, T(0.0))
    end
end

# ===== Constructors

function PSO(prob::Problem{fType,S}; numParticles = 100, inertiaRange = (0.1, 1.1), minNeighborFrac = 0.25, 
    selfAdjustWeight = 1.49, socialAdjustWeight = 1.49, initMethod::Symbol = :Uniform,
    updateMethod::Symbol = :MATLAB) where {S,fType}

    # Error checking
    length(inertiaRange) == 2 || throw(ArgumentError("inertiaRange must be of length 2."))
    minNeighborFrac > 0       || throw(ArgumentError("minNeighborFrac must be > 0."))

    # Type info
    T = typeof(1.0)
    nIRange = (T(inertiaRange[1]), T(inertiaRange[2]))

    # Call constructor
    return PSO{T}(prob, numParticles, initMethod, updateMethod, nIRange, T(minNeighborFrac),
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

    # Set optimizer state 
    pso.state    = 1

    # Call callback function 
    if opts.callback !== nothing
        opts.callback(pso, opts)
    end

    # Print Status
    if opts.display
        printStatus(pso.swarm, 0.0, 0, 0)
    end

    return nothing
end

function iterate!(pso::PSO, opts::Options)

    # Initialize time and iteration counter
    pso.t0  		= time()
    pso.stallT0 	= pso.t0
    pso.fStall 		= Inf

    # Compute minimum neighborhood size
    minNeighborSize = max(2, floor(length(pso.swarm)*pso.minNeighborFrac))

    # Begin loop
    exitFlag = 0
    while exitFlag == 0

        # Update iteration counter 
        pso.iters += 1

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
        if pso.fStall - pso.swarm.b > opts.funcTol
            pso.fStall 		= pso.swarm.b 
            pso.stallIters 	= 0
            pso.stallT0 	= time()
        else
            pso.stallIters 	+= 1
        end

        # Stopping criteria
        if pso.stallIters >= opts.maxStallIters
            exitFlag 	= 1
	    pso.state 	= 3
        elseif pso.iters >= opts.maxIters
            exitFlag 	= 2
	    pso.state 	= 3
        elseif pso.swarm.b <= opts.objLimit
            exitFlag 	= 3
	    pso.state 	= 3
        elseif time() - pso.stallT0 >= opts.maxStallTime 
            exitFlag 	= 4
	    pso.state 	= 3
        elseif time() - pso.t0 >= opts.maxTime 
            exitFlag 	= 5
	    pso.state 	= 3
    	else
            pso.state 	= 2
        end

        # Output Status
        if opts.display && pso.iters % opts.displayInterval == 0
            printStatus(pso.swarm, time() - pso.t0, pso.iters, pso.stallIters)
        end

        # Call callback function
        if opts.callback !== nothing
            opts.callback(pso, opts)
        end
    end

    # Set results
    setResults!(pso, exitFlag)

    return nothing
end

function setResults!(pso::PSO, exitFlag::Int)
    pso.results.fbest = pso.swarm.b 
    pso.results.xbest .= pso.swarm.d
    pso.results.iters = pso.iters 
    pso.results.time = time() - pso.t0
    pso.results.exitFlag = exitFlag

    return nothing
end
