
struct MS_PSO{T,S,fType} <: Optimizer 

    # Optimization problem 
    prob::Problem{fType,S}

    # Vector of particle swarms
    swarmVec::Vector{Swarm{T}}

    # Settings 
    initMethod::Symbol 
    updateMethod::Symbol

    # PSO specific parameters 
    inertiaRange::Tuple{T,T}
    minNeighborFrac::T
    selfAdjustWeight::T
    socialAdjustWeight::T

    # Results data 
    results::Results{T}

    function MS_PSO{T}(prob::Problem{fType,S}, numParticlesPerSwarm::Integer, 
        numSwarms::Integer, initMethod::Symbol, updateMethod::Symbol, 
        inertiaRange::Tuple{T,T}, minNeighborFrac::T, selfAdjustWeight::T, 
        socialAdjustWeight::T) where {T,S,fType}

        # Initialize Swarm vector
        N = length(prob.LB)
        swarmVec = [Swarm{T}(N, numParticlesPerSwarm) for i in 1:numSwarms]

        # Initialize results
        results = Results{T}(undef, N)

        # Check that init and update methods exist
        if initMethod != :Uniform && initMethod != :LogisticsMap
            throw(ArgumentError("$initMethod is not a PSO initialization method."))
        end
        if updateMethod != :MATLAB 
            throw(ArgumentError("$updateMethod is not a PSO update method."))
        end

        return new{T,S,fType}(prob, swarmVec, initMethod, updateMethod, 
                        inertiaRange, minNeighborFrac, selfAdjustWeight, 
                        socialAdjustWeight, results)
    end
end

# ===== Constructors

function MS_PSO(prob::Problem{fType,S}; numParticlesPerSwarm = 100, numSwarms = 4, 
    inertiaRange = (0.1, 1.1), minNeighborFrac = 0.25, selfAdjustWeight = 1.49, 
    socialAdjustWeight = 1.49, initMethod::Symbol = :Uniform, 
    updateMethod::Symbol = :MATLAB) where {S,fType}

    # Error checking
    length(inertiaRange) == 2 || throw(ArgumentError("inertiaRange must be of length 2."))
    minNeighborFrac > 0       || throw(ArgumentError("minNeighborFrac must be > 0."))

    # Type info
    T = typeof(1.0)
    nIRange = (T(inertiaRange[1]), T(inertiaRange[2]))

    # Call constructor
    return MS_PSO{T}(prob, numParticlesPerSwarm, numSwarms, initMethod, updateMethod, 
        nIRange, T(minNeighborFrac), T(selfAdjustWeight), T(socialAdjustWeight))
end

# ===== Methods

function _optimize!(mspso::MS_PSO, opts::Options)

    # Initialize PSO
    initialize!(mspso, opts)

    # Iterate PSO
    iterate!(mspso, opts)

    # Display results and return 
    if opts.display
        display(mspso.results)
    end

    return mspso.results
end

function initialize!(mspso::MS_PSO, opts::Options)
    @inbounds begin
    # Initialize position and velocities
    if mspso.initMethod == :Uniform
        for i in 1:length(mspso.swarmVec)
            uniformInitialization!(mspso.swarmVec[i], mspso.prob, opts)
        end
    elseif mspso.initMethod == :LogisticsMap
        for i in 1:length(mspso.swarmVec)
            logisticsMapInitialization!(mspso.swarmVec[i], mspso.prob, opts)
        end
    else
        throw(ArgumentError("PSO initialization method is not implemented."))
    end

    # Iterate through each swarm
    # Note: This implementation is not entirely optimal as we will 
    # be waiting for each swarm to finish initializing before starting the next
    for i in 1:length(mspso.swarmVec)
        # Evaluate Objective Function
        feval!(mspso.swarmVec[i], mspso.prob.f, opts; init = true)

        # Set swarm global best obj. func. value to Inf 
        mspso.swarmVec[i].b = Inf

        # Set global best 
        setGlobalBest!(mspso.swarmVec[i])

        # Initialize neighborhood size
        mspso.swarmVec[i].n = max(2, floor(length(mspso.swarmVec[i])*mspso.minNeighborFrac))

        # Initial inertia
        if mspso.inertiaRange[2] > 0
            mspso.swarmVec[i].w = (mspso.inertiaRange[2] > mspso.inertiaRange[1]) ?
                                mspso.inertiaRange[2] : mspso.inertiaRange[1]
        else
            mspso.swarmVec[i].w = (mspso.inertiaRange[2] < mspso.inertiaRange[1]) ?
                                mspso.inertiaRange[2] : mspso.inertiaRange[1]
        end

        # Initialize self and social adjustment 
        mspso.swarmVec[i].y₁ = mspso.selfAdjustWeight
        mspso.swarmVec[i].y₂ = mspso.socialAdjustWeight
    end

    # Print Status # NEED TO IMPLEMENT A NEW printStatus for Vector{Swarm}
    if opts.display
        printStatus(mspso.swarmVec, 0.0, 0, 0)
    end

    end
    return nothing
end

function reset!(mspso::MS_PSO, opts::Options, resetIdx)
    # GH: Only allow a swarm to reset if another swarm with a better solution exists
    # If we don't do this, we can end up throwing away our best found solution.
    if mspso.swarmVec[resetIdx].b > getBestF(mspso)
        if mspso.initMethod == :Uniform
            uniformInitialization!(mspso.swarmVec[resetIdx], mspso.prob, opts)
        elseif mspso.initMethod == :LogisticsMap
            logisticsMapInitialization!(mspso.swarmVec[resetIdx], mspso.prob, opts)
        else
            throw(ArgumentError("PSO initialization method is not implemented."))
        end
        # Evaluate Objective Function
        feval!(mspso.swarmVec[resetIdx], mspso.prob.f, opts; init = true)

        # Set swarm global best obj. func. value to Inf 
        mspso.swarmVec[resetIdx].b = Inf

        # Set global best 
        setGlobalBest!(mspso.swarmVec[resetIdx])

        # Initialize neighborhood size
        mspso.swarmVec[resetIdx].n = max(2, floor(length(mspso.swarmVec[resetIdx])*mspso.minNeighborFrac))
    end
    return nothing
end

function iterate!(mspso::MS_PSO, opts::Options)

    # Initialize time and iteration counter
    t0  = time()
    stallT0 = t0*ones(length(mspso.swarmVec))
    iters = 0
    stallIters = fill!(Vector{Int}(undef, length(mspso.swarmVec)), 0)
    hasStalled = BitVector(undef, length(mspso.swarmVec))
    fStall = Inf*ones(length(mspso.swarmVec))

    # Compute minimum neighborhood size
    minNeighborSize = max(2, floor(length(mspso.swarmVec[1])*mspso.minNeighborFrac))

    # Begin loop
    exitFlag = 0
    while exitFlag == 0

        # Update iteration counter 
        iters += 1

        # Prepare all swarms for evaluation
        for i in 1:length(mspso.swarmVec)
            # Update swarm velocities 
            updateVelocities!(mspso.swarmVec[i])

            # Update positions 
            step!(mspso.swarmVec[i])

            # Enforce Bounds
            if !(all(isinf.(mspso.prob.LB)) && all(isinf.(mspso.prob.UB)))
                enforceBounds!(mspso.swarmVec[i], mspso.prob.LB, mspso.prob.UB)
            end
        end

        # Evaluate objective function for all swarms 
        # Will need to experiment with multithreadding here!
        for i in 1:length(mspso.swarmVec)
            feval!(mspso.swarmVec[i], mspso.prob.f, opts)
        end

        # Update all swarms
        for i in 1:length(mspso.swarmVec)
            # Update global best
            flag::Bool = setGlobalBest!(mspso.swarmVec[i])

            # Update stall counter and neighborhood
            if flag
                mspso.swarmVec[i].c = max(0, mspso.swarmVec[i].c - 1)
                mspso.swarmVec[i].n = minNeighborSize
            else 
                mspso.swarmVec[i].c += 1
                mspso.swarmVec[i].n = min(mspso.swarmVec[i].n + minNeighborSize, 
                                      length(mspso.swarmVec[i]) - 1)
            end

            # Update inertia
            if mspso.swarmVec[i].c < 2
                mspso.swarmVec[i].w *= 2.0
            elseif mspso.swarmVec[i].c > 5
                mspso.swarmVec[i].w /= 2.0
            end
            
            # Ensure inertia is in inertia range 
            if  mspso.swarmVec[i].w > mspso.inertiaRange[2]
                mspso.swarmVec[i].w = mspso.inertiaRange[2]
            elseif mspso.swarmVec[i].w < mspso.inertiaRange[1]
                mspso.swarmVec[i].w = mspso.inertiaRange[1]
            end 

            #Check if two swarms are too close
            distancecheck!(mspso, opts)

            # Track stalling
            if fStall[i] - mspso.swarmVec[i].b > opts.funcTol
                fStall[i] = mspso.swarmVec[i].b 
                stallIters[i] = 0
                stallT0[i] = time()
            else
                stallIters[i] += 1
            end

            # Check for stalling 
            if stallIters[i] >= opts.maxStallIters
                hasStalled[i] = true
                interMingle!(mspso, opts, i) # GH: Updated arguments pased to interMingle! 
            end

            
       end

        # Stopping criteria
        if all(hasStalled)
            exitFlag = 1
        elseif iters >= opts.maxIters
            exitFlag = 2
        elseif any([mspso.swarmVec[i].b <= opts.objLimit for i in 1:length(mspso.swarmVec)])
            exitFlag = 3
        elseif all([time() - stallT0[i] >= opts.maxStallTime for i in 1:length(mspso.swarmVec)])
            exitFlag = 4
        elseif time() - t0 >= opts.maxTime 
            exitFlag = 5
        end

        # Output Status
        if opts.display && iters % opts.displayInterval == 0
            printStatus(mspso.swarmVec, time() - t0, iters, stallIters)
        end
    end

    # Set results
    setResults!(mspso, iters, time() - t0, exitFlag)

    return nothing
end

# GH: Not quite sure I understand where you ment for this block of code to be.
# As is, this block is not within any function, so swarm, mspso, opts, and reset bestIdx
# are not defined. Commenting out block for now. 

# Just for you to get a better understanding of how Julia works, when code like this is written
# outside of a function, it's in the "global scope" and only has access to variables defined in 
# the global scope (variables in the global scope are usually defined in the main script you write
# or are functions exported from the packages you import). 
# Likewise, when code is in a function, it only has access to the varables within that functions scope
# (this is why I needed to modify the interMingle! function, to ass mspso and opts to the functions 
# scope). In Julia, in general, we want to put as much as possible within a function, because
# only code within functions is compiled (and therefor runs fast). Code in the global scope
# is NOT compiled. 

#Calculate the distances between best for each swarm and reset if too close
function distancecheck!(mspso::MS_PSO, opts::Options)
    #initialize distance check matrix
    r=zeros(Float64,4,4)
@inbounds begin
for i=1:4, j=1:4 
    DistVec=(mspso.swarmVec[j].d-mspso.swarmVec[i].d)
    r[i,j]=sqrt(DistVec[1]^2+DistVec[2]^2)
     if r[i,j]==0.0
        r[i,j]=10.0
     elseif r[i,j]<2.0
        resetIdx=i
     reset!(mspso,opts,resetIdx)
    end
end
end
end


# Could intermingle stalled swarm with best swarm or random swarm. Choosing random swarm for now!
function interMingle!(mspso::MS_PSO, opts::Options, stallIdx) where {T} # GH: Updated intermingle arguments to include MS_PSO object instead of just the swarm vector and options
    # resetting stalled swarm
    # GH: Removed extra resetIdx variable for simplicity.
    reset!(mspso, opts, stallIdx) # GH: Resetting each swarm every time max stall iters is hit results in failiur to converge... We'll need to update how this is done
                                  # If you comment out the call to reset! here, the tests should pass.
    return nothing
end

function printStatus(swarmVec::Vector{Swarm{T}}, time, iter, stallCount::Vector{Int}) where {T}
    bestF = minimum([swarmVec[i].b for i in 1:length(swarmVec)])
    fspec1 = FormatExpr("Time Elapsed: {1:f} sec, Iteration Number: {2:d}, Function Evaluations: {3:d}")
    fspec2 = FormatExpr("Max Stall Iterations: {1:d}, Global Best: {2:e}")
    printfmtln(fspec1, time, iter, (iter + 1)*length(swarmVec[1])*length(swarmVec))
    printfmtln(fspec2, maximum([stallCount[i] for i in 1:length(swarmVec)]), bestF)
    println(" ")
end

function printStatus(swarmVec::Vector{Swarm{T}}, time, iter, stallCount::Int) where {T}
    bestF = minimum([swarmVec[i].b for i in 1:length(swarmVec)])
    fspec1 = FormatExpr("Time Elapsed: {1:f} sec, Iteration Number: {2:d}, Function Evaluations: {3:d}")
    fspec2 = FormatExpr("Max Stall Iterations: {1:d}, Global Best: {2:e}")
    printfmtln(fspec1, time, iter, (iter + 1)*length(swarmVec[1])*length(swarmVec))
    printfmtln(fspec2, stallCount, bestF)
    println(" ")
end

# Function to get the best objective function value from all swarms
function getBestF(mspso::MS_PSO)
    bestF = Inf 
    for i in 1:length(mspso.swarmVec)
        if mspso.swarmVec[i].b < bestF 
            bestF = mspso.swarmVec[i].b
        end
    end
    return bestF
end

function setResults!(mspso::MS_PSO, iters::Int, Δt::AbstractFloat, exitFlag::Int)
    bestIdx = -1
    bestF = Inf
    for i in 1:length(mspso.swarmVec)
        if mspso.swarmVec[i].b < bestF 
            bestF = mspso.swarmVec[i].b
            bestIdx = i
        end
    end
    mspso.results.fbest = bestF
    mspso.results.xbest .= mspso.swarmVec[bestIdx].d
    mspso.results.iters = iters 
    mspso.results.time = Δt 
    mspso.results.exitFlag = exitFlag

    return nothing
end