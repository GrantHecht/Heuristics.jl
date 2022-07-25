using Heuristics, StaticArrays, BenchmarkTools

include("c:\\Users\\Owner\\.julia\\dev\\Heuristics\\test\\testProblems.jl")

# Size of problem 
N = 8
M = 500

# Define problem 
LB = -10 .* ones(N)
UB =  10 .* ones(N)
prob = Problem(rastriginfunc, LB, UB)
opts = Options(display=true, maxIters=1000, maxStallIters = 1000, funcTol = 1e-8)

# Instantiate MPSO
mspso = MS_PSO(prob; numParticlesPerSwarm = M)

# optimize
res = optimize!(mspso, opts)