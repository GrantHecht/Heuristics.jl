using Heuristics, StaticArrays, BenchmarkTools

include("./test/testProblems.jl")

# Size of problem 
N = 4
M = 100

# Define problem 
LB = -15 .* ones(N)
UB =  30 .* ones(N)
prob = Problem(ackleyfunc, LB, UB)
opts = Options(display=false, maxIters=1000, maxStallIters = 1000, funcTol = 1e-8)

# Instantiate MPSO
mpso = MPSO(prob; numParticles = M)

# optimize
res = optimize!(mpso, opts)