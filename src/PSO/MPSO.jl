# Modified PSO

mutable struct MPSO{T,S,fType} <: Optimizer

    # Optimization problem
    prob::Problem{fType,S}

    # Swarm of particles
    swarm::Swarm{T}

    # Settings
    minNeighborFrac::T

    # MPSO specific parameters 
    MFDstar::T
    α::T
    c1::T 
    c2::T
    sc::Int
    fc::Int

    g::T # Max value dilation parameter can take
    ψ::T # Shape parameter
    pm::T

    a::T
    MFD::T
    fAvg::T
    ρ::T

    # Results data
    results::Results{T}

    function MPSO{T}(prob::Problem{fType,S}, numParticles::Integer, α, c1, c2, sc, fc, g, ψ, pm) where {T,S,fType}

        # Initialize Swarm 
        N = length(prob.LB)
        swarm = Swarm{T}(N, numParticles)

        # Compute and set pm if necessary
        if pm == -1.0
            n = length(prob.LB)
            if n > 1000
                pm = 0.025
            elseif n > 101
                pm = 0.075
            elseif n > 11
                pm = 0.15
            elseif n > 5
                pm = 0.35
            else
                pm = 0.5
            end
        end

        # Initialize results
        results = Results{T}(undef, N)

        return new{T,S,fType}(prob, swarm, T(0.25), T(2.8e-6), T(α), T(c1), T(c2), sc, fc, T(g), T(ψ), T(pm), T(NaN), T(NaN), T(NaN), T(1.0), results)
    end
end

# ===== Constructors

function MPSO(prob::Problem{fType,S}; 
    numParticles = 100, α = 0.2, c1 = 1.49, c2 = 1.49, sc = 15, fc = 5, g = 10000.0, ψ = 1.0, pm = -1.0) where {S,fType}

    # Type info
    T = typeof(1.0)

    # Call constructor
    return MPSO{T}(prob, numParticles, α, c1, c2, sc, fc, g, ψ, pm)
end

# ===== Methods 

function _optimize!(mpso::MPSO, opts::Options)
    # Initialize MPSO
    initialize!(mpso, opts)

    # Iterate MPSO
    iterate!(mpso, opts)

    # Display results and return 
    if opts.display
        display(mpso.results)
    end

    return mpso.results
end

function initialize!(mpso::MPSO, opts::Options)
    # Initialize with logistics map
    logisticsMapInitialization!(mpso.swarm, mpso.prob, opts)

    # Evaluate objective function 
    feval!(mpso.swarm, mpso.prob.f, opts; init = true)

    # Set global best
    mpso.swarm.b = Inf
    setGlobalBest!(mpso.swarm)

    # Set inertia 
    mpso.swarm.w = 0.9

    # Print Status
    if display
        #printStatus(mpso.swarm, 0.0, 0.0)
    end

    return nothing
end

function iterate!(mpso::MPSO, opts::Options)

    # Initialize time and iteration counter
    t0 = time()
    iters = 0
    stallIters = 0
    fStall = Inf

    # Begin loop 
    exitFlag    = 0
    prevSucc    = true
    nSeq        = 0
    while exitFlag == 0
        # Compute average fitness and MFD
        mpso.fAvg = computeAverangeFitness(mpso.swarm); 
        mpso.MFD  = computeMFD(mpso.swarm)

        # Update position of particles based on MFD
        bestIdx = updateGlobalBestParticlePosition!(mpso)
        if mpso.MFD > mpso.MFDstar || iters == 0
            updateParticlePositions!(mpso, bestIdx)
        else
            throw(ErrorException("Fuck"))
            waveletMutateOrReinit!(mpso, bestIdx, iters, opts.maxIters)
        end

        # Enfore Bounds 
        if !(all(isinf.(mpso.prob.LB)) && all(isinf.(mpso.prob.UB)))
            enforceBounds!(mpso.swarm, mpso.prob.LB, mpso.prob.UB)
        end

        # Evaluate objective function
        feval!(mpso.swarm, mpso.prob.f, opts)

        # Update global best 
        succ::Bool = setGlobalBest!(mpso.swarm)

        # Update ρ
        if prevSucc == succ 
            nSeq += 1
            if succ == true && nSeq > mpso.sc 
                mpso.ρ *= 2.0
            elseif succ == false && nSeq > mpso.fc 
                mpso.ρ *= 0.5
            end
        else
            prevSucc = succ 
            nSeq = 0
        end

        # Update ω
        mpso.swarm.w = computeInertiaWeight(iters, opts.maxIters, mpso.α)

        # Increment iteration counter 
        iters += 1

        # Stall tracking
        if fStall - mpso.swarm.b > opts.funcTol 
            fStall = mpso.swarm.b
            stallIters = 0
        else
            stallIters += 1
        end

        # Stopping criteria 
        if stallIters >= opts.maxStallIters
            exitFlag = 1
        elseif iters >= opts.maxIters 
            exitFlag = 2
        elseif mpso.swarm.b <= opts.objLimit
            exitFlag = 3
        end

        # Output status 
        if opts.display && iters % opts.displayInterval == 0
            printStatus(mpso.swarm, time() - t0, iters, mpso.MFD, mpso.fAvg)
        end
    end

    # Set results
    setResults!(mpso, iters, time() - t0, exitFlag)

    return nothing
end

function setResults!(mpso::MPSO, iters::Int, Δt::AbstractFloat, exitFlag::Int)
    mpso.results.fbest  = mpso.swarm.b
    mpso.results.xbest .= mpso.swarm.d 
    mpso.results.iters  = iters
    mpso.results.time   = Δt 
    mpso.results.exitFlag = exitFlag

    return nothing
end

function computeInertiaWeight(iters, maxIters, α)
    ω = 0.9
    if iters > α*maxIters 
        ω = 1.0 / (1.0 + exp((10.0*iters - 2.0*maxIters) / maxIters)) + 0.4
    end
    return ω
end

function updateGlobalBestParticlePosition!(mpso::MPSO)
    @inbounds begin
        # Find best particle
        f = Inf
        bestIdx = -1 
        @simd for i in 1:length(mpso.swarm)
            if mpso.swarm[i].fp < f
                f = mpso.swarm[i].fp 
                bestIdx = i
            end 
        end

        # Update position 
        @simd for i in 1:length(mpso.swarm.d)
            mpso.swarm[bestIdx].v[i] = -mpso.swarm[bestIdx].x[i] + mpso.swarm.d[i] + mpso.swarm.w*mpso.swarm[bestIdx].v[i] + 
                                    mpso.ρ*(1.0 - 2.0*rand())
            mpso.swarm[bestIdx].x[i] += mpso.swarm[bestIdx].v[i]
        end
    end
    return bestIdx
end

function updateParticlePositions!(mpso::MPSO, bestIdx)
    @inbounds begin
        @simd for i in 1:length(mpso.swarm)
            if i != bestIdx
                @simd for j in 1:length(mpso.swarm.d)
                    mpso.swarm[i].v[j] = mpso.swarm.w*mpso.swarm[i].v[j] + 
                                         mpso.c1*rand()*(mpso.swarm[i].p[j] - mpso.swarm[i].x[j]) +
                                         mpso.c2*rand()*(mpso.swarm.d[j] - mpso.swarm[i].x[j])
                    mpso.swarm[i].x[j] += mpso.swarm[i].v[j]
                end
            end
        end
    end
end

function waveletMutateOrReinit!(mpso, bestIdx, iters, maxIters)
    @inbounds begin
        updateDilationParam!(mpso, iters, maxIters)
        @simd for i in 1:length(mpso.swarm)
            if i != bestIdx
                if mpso.swarm[i].fx <= mpso.fAvg # Mutate!
                    @simd for j in 1:length(mpso.swarm.d)
                        if rand() <= mpso.pm
                            σ = computeSigma(mpso)
                            if σ > 0
                                mpso.swarm[i].x[j] += σ*(mpso.prob.UB[j] - mpso.swarm[i].x[j])
                            else
                                mpso.swarm[i].x[j] += σ*(mpso.swarm[i].x[j] - mpso.prob.LB[j])
                            end
                        end
                    end
                else
                    reInitializeParticle!(mpso::MPSO, i)
                end
            end
        end
    end
end

function updateDilationParam!(mpso::MPSO, iters, maxIters)
    mpso.a = exp(-log(mpso.g)*(1.0 - iters / maxIters)^mpso.ψ) + log(mpso.g)
    return nothing
end

function computeSigma(mpso::MPSO)
    ϕ = mpso.c1 + mpso.c2
    return cos(5.0 * ϕ / mpso.a)*exp(-(ϕ / mpso.a)^2 / 2.0) / sqrt(mpso.a)
end

function reInitializeParticle!(mpso::MPSO, idx)
    N   = length(mpso.swarm.d)
    LB  = mpso.prob.LB
    UB  = mpso.prob.UB

    fixedPointTol = 1e-14
    maxPert = 1e-12
    lMapIters = 3000
    @inbounds begin
        @simd for j in 1:N
            mpso.swarm[idx].x[j] = 0.4567 + 2.0*(rand() - 0.5)*maxPert
        end
        @simd for j in 1:N
            @simd for k in 1:lMapIters
                val = mpso.swarm[idx].x[j]
                if val < fixedPointTol || 
                    abs(val - 0.25) < fixedPointTol || 
                    abs(val - 0.50) < fixedPointTol ||
                    abs(val - 0.75) < fixedPointTol ||
                    abs(val - 1.00) < fixedPointTol

                    mpso.swarm[idx].x[j] += maxPert*rand()
                end
                mpso.swarm[idx].x[j] = lMap(mpso.swarm[idx].x[j])
                if isinf(mpso.swarm[idx].x[j])
                    throw(ErrorException("Inf or Nan encountered during logistic map reinitialization."))
                end
            end
        end
    end

    # Scale particle positions 
    @simd for j in 1:N
        mpso.swarm[idx].x[j] = LB[j] + (UB[j] - LB[j])*mpso.swarm[idx].x[j]
    end
end
