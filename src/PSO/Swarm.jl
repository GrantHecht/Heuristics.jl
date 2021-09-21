mutable struct Swarm{T<:AbstractFloat}

    # Vector of particles
    particles::Vector{Particle{T}}

    # Preallocated vector of UInt16 for neighborhood selection
    nVec::Vector{Int}

    # Global best objective function value
    b::T

    # Location of global best objective function value 
    d::Vector{T}

    n::Int      # Neighborhood size 
    w::T        # Inertia
    c::Int      # Adaptive inertia counter

    y₁::T   # Self adjustment weight
    y₂::T   # Social adjustment weight

end

function Swarm{T}(::UndefInitializer) where {T}
        return Swarm{T}(Vector{Particle{T}}(undef, 0),
            Vector{Int}(undef, 0), T(0.0), Vector{T}(undef, 0), 
            Int(0), T(0.0), Int(0), T(0.0), T(0.0))
end

function Swarm{T}(nDims::Integer, nParticles::Integer) where {T}
    if nDims < 0
        throw(ArgumentError("nDims must be greater than 0.")) 
    end
    if nParticles < 0
        throw(ArgumentError("nParticles must be greater than 0."))
    end
    return Swarm{T}(Vector{Particle{T}}([Particle{T}(nDims) for i in 1:nParticles]),
        Vector{Int}(1:nParticles), T(0.0), Vector{T}(undef, nDims), 
        Int(0), T(0.0), Int(0), T(0.0), T(0.0))
end


# ===== Methods

Base.length(s::Swarm) = length(s.particles)

function Base.getindex(s::Swarm, i::Int)
    1 <= i <= length(s.particles) || throw(BoundsError(s, i))
    return s.particles[i]
end

function Base.setindex!(s::Swarm, v, i::Int)
    1 <= i <= length(s.particles) || throw(BoundsError(s, i))
    s.particles[i] = v
end

function feval!(s::Swarm, f, opts::Options; init = false)
    # Evaluate objective functions
    @inbounds begin
        if opts.useParallel
            @qthreads for i in 1:length(s)
                s[i].fx = f(s[i].x)
            end
        else
            @simd for i in 1:length(s) 
                s[i].fx = f(s[i].x)
            end
        end
    end

    # Check objective function values if desired 
    if opts.funValCheck
        fValCheck(s)
    end

    # Update each particles best objective function value and its location
    if !init
        @inbounds @simd for i in 1:length(s)
            if s[i].fx < s[i].fp
                s[i].p .= s[i].x
                s[i].fp = s[i].fx
            end
        end
    else
        @inbounds @simd for i in 1:length(s)
            s[i].fp = s[i].fx
        end
    end
end

function setGlobalBest!(s::Swarm)
    updated = false
    @inbounds begin
        @simd for i in 1:length(s)
            if s[i].fp < s.b 
                s.b = s[i].fp
                s.d .= s[i].p
                updated = true
            end
        end
    end
    return updated
end

function fValCheck(s::Swarm)
    @inbounds begin
        @simd for i in 1:length(s)
            if isinf(s[i].fx) || isnan(s[i].fx)
                throw(ErrorException("Objective function returned Inf or NaN!"))
            end
        end      
    end
end

function updateVelocities!(s::Swarm)
    n = length(s.d)
    m = length(s)
    @inbounds begin
        @simd for i in 1:m

            # Shuffle vector containing integers 1:n
            # first m != i will be neighborhood
            shuffle!(s.nVec)

            # Determine fbest(S)
            fbest = Inf
            best = 0
            incr = 0
            for j in 1:s.n
                s.nVec[j] == i ? incr = 1 : () 
                k = j + incr
                if s[k].fp < fbest
                    fbest = s[k].fp
                    best = k
                end
            end

            # Update i's velocity 
            @simd for j in 1:n
                s[i].v[j] = s.w*s[i].v[j] + s.y₁*rand()*(s[i].p[j] - s[i].x[j]) + 
                                s.y₂*rand()*(s[best].p[j] - s[i].x[j])
            end 
        end
    end
end

function step!(s::Swarm)
    @inbounds @simd for i in 1:length(s)
        s[i].x .= s[i].x .+ s[i].v 
    end
end

function enforceBounds!(s::Swarm, LB, UB)
    @inbounds begin
        @simd for i in 1:length(s)
            @simd for j in 1:length(s.d)
                if s[i].x[j] > UB[j]
                    s[i].x[j] = UB[j]
                    s[i].v[j] = 0.0
                elseif s[i].x[j] < LB[j]
                    s[i].x[j] = LB[j]
                    s[i].v[j] = 0.0
                end
            end
        end
    end
end

function computeAverangeFitness(s::Swarm)
    @inbounds begin
        sum = 0.0
        @simd for i in 1:length(s)
            sum += s[i].fx 
        end
    end
    return sum / length(s)
end

function computeMFD(s::Swarm)
    @inbounds begin
        MFD = -Inf
        @simd for i in 1:length(s)
            sum = 0.0
            @simd for j in 1:length(s.d)
                #sum += (s[i].p[j] - s[i].x[j])^2
                sum += (s.d[j] - s[i].x[j])^2
            end
            sqrtTerm = sqrt(sum / length(s.d))
            if sqrtTerm > MFD 
                MFD = sqrtTerm
            end
        end
    end
    return MFD
end

function printStatus(s::Swarm, time::AbstractFloat, iter::Int, stallCount::Int)
    fspec1 = FormatExpr("Time Elapsed: {1:f} sec, Iteration Number: {2:d}, Function Evaluations: {3:d}")
    fspec2 = FormatExpr("Stall Iterations: {1:d}, Global Best: {2:e}")
    printfmtln(fspec1, time, iter, (iter + 1)*length(s))
    printfmtln(fspec2, stallCount, s.b)
    println(" ")
end

function printStatus(s::Swarm, time::AbstractFloat, iter::Int, MFD::AbstractFloat, fAvg::AbstractFloat)
    fspec1 = FormatExpr("Time Elapsed: {1:f} sec, Iteration Number: {2:d}, Function Evaluations: {3:d}")
    fspec2 = FormatExpr("MFD: {1:e}, Avg. Fitness: {2:e}, Global Best: {3:e}")
    printfmtln(fspec1, time, iter, (iter + 1)*length(s))
    printfmtln(fspec2, MFD, fAvg, s.b)
    println(" ")
end

# Initializes position and velocities of particles by sampling from a 
# uniform distribution
function uniformInitialization!(swarm::Swarm, prob::Problem, opts::Options)

    # Get N: Number of diamensions and M: Swarm Size 
    M = length(swarm)
    N = length(swarm[1])

    # Get Boundary Constraints
    LB  = prob.LB
    UB  = prob.UB

    # Check if initial bounds on positions have been set
    useInitBnds = false
    if length(opts.iLB) == N && length(opts.iUB) == N
        useInitBnds = true
        iLB = opts.iLB            
        iUB = opts.iUB
    end

    # Initialize particle positions and velocities
    @inbounds begin
        @simd for d in 1:N
            # Get local bounds for d-axis
            lLB = useInitBnds ? (LB[d] < iLB[d] ? iLB[d] : LB[d]) : LB[d]
            lUB = useInitBnds ? (UB[d] > iUB[d] ? iUB[d] : UB[d]) : UB[d]
            @simd for p in 1:M
                # Position information
                swarm[p].x[d] = lLB + (lUB - lLB)*rand()
                swarm[p].p[d] = swarm[p].x[d]

                # Velocity 
                r = useInitBnds ? min(lUB-lLB,UB[d]-LB[d]) : lUB - lLB 
                swarm[p].v[d] = -r + 2*r*rand()
            end
        end
    end
    
    return nothing
end

# Initializes position and velocities of particles using logistic map
function logisticsMapInitialization!(swarm::Swarm, prob::Problem, opts::Options)

    # Get N: Number of diamensions and M: Swarm Size 
    M = length(swarm)
    N = length(swarm[1])

    # Get Boundary Constraints
    LB  = prob.LB
    UB  = prob.UB

    # Check if initial bounds on positions have been set
    useInitBnds = false
    if length(opts.iLB) == N && length(opts.iUB) == N
        useInitBnds = true
        iLB = opts.iLB            
        iUB = opts.iUB
    end

    # Logistics map initialization
    fixedPointTol = 1e-14
    maxPert = 1e-12
    lMapIters = 3000
    @inbounds begin
        @simd for j in 1:M
            @simd for k in 1:N
                swarm[j].x[k] = 0.4567 + 2*(rand() - 0.5)*maxPert
            end
        end
        @simd for i in 1:lMapIters
            @simd for j in 1:M 
                @simd for k in 1:N
                    val = swarm[j].x[k]
                    if val < fixedPointTol || 
                       abs(val - 0.25) < fixedPointTol || 
                       abs(val - 0.50) < fixedPointTol ||
                       abs(val - 0.75) < fixedPointTol ||
                       abs(val - 1.00) < fixedPointTol

                       swarm[j].x[k] += maxPert*rand()
                    end
                    swarm[j].x[k] = lMap(swarm[j].x[k])
                    if isinf(swarm[j].x[k]) 
                        throw(ErrorException("Inf or NaN"))
                    end
                end
            end
        end
    end

    # Scale particle positions and initialize velocities
    @inbounds begin
        @simd for d in 1:N
            # Get local bounds for d-axis
            lLB = useInitBnds ? (LB[d] < iLB[d] ? iLB[d] : LB[d]) : LB[d]
            lUB = useInitBnds ? (UB[d] > iUB[d] ? iUB[d] : UB[d]) : UB[d]
            @simd for p in 1:M
                # Position information
                swarm[p].x[d] = lLB + (lUB - lLB)*swarm[p].x[d]
                swarm[p].p[d] = swarm[p].x[d]

                # Velocity 
                r = useInitBnds ? min(lUB-lLB,UB[d]-LB[d]) : lUB - lLB 
                swarm[p].v[d] = -r + 2*r*rand()
            end
        end
    end
    
    return nothing
end