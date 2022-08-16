
using Heuristics, StaticArrays, BenchmarkTools

# Objective function
include("./../test/testProblems.jl")

function main()

    # Size of problem
    N = 8 # Length of the decision vector
    M = 500 # Number of particles

    # Define problem: N is number of dims 
    LB = -10 .* ones(N)
    UB =  10 .* ones(N)
    prob = Problem(rastriginfunc, LB, UB)
    opts = Options(display=true)

    # Initialize Optimizer: M is swarm size!
    optimizer = PSO(prob; numParticles = M)

    # solve 
    #Heuristics.initialize!(optimizer, opts)
    #@btime Heuristics.iterate!($optimizer, $opts)
    optimize!(optimizer, opts)

    return nothing
end

main()