using Heuristics, StaticArrays, BenchmarkTools

include("./../test/testProblems.jl")

# Size of problem 
N = 8 # Length of the decision vector
M = 500 # Number of particles

# Define problem 
LB = -10 .* ones(N)
UB =  10 .* ones(N)
prob = Problem(rastriginfunc, LB, UB)
#opts = Options(;display=true, maxIters=1000, maxStallIters = 1000, funcTol = 1e-8)
opts = Options(;display = true, maxIters=1000, maxStallIters = 1000)

# Instantiate MPSO
mspso = MS_PSO(prob; numParticlesPerSwarm = M)

# optimize
res = optimize!(mspso, opts)