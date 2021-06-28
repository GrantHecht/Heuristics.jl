
using Heuristics, StaticArrays

# Objective function
f(x) = sum(x.^2)

function main()

    # Size of problem
    N = 5
    M = 50

    # Define problem: N is number of dims 
    LB = -5 .* ones(N)
    UB =  5 .* ones(N)
    prob = Problem{N}(f, LB, UB)

    # Initialize Optimizer: M is swarm size!
    optimizer = PSO{M}(prob;)

    # solve 
    optimize!(optimizer)

    return nothing
end

main()