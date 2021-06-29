using Heuristics

include("testProblems.jl")

# Sphere Function
N = 14
M = 200
LB = -5 .* ones(N)
UB = 5 .* ones(N)
prob = Problem{N}(sphereFunc, LB, UB)
opts = Options(;display = false)
pso = PSO{M}(prob)
res = optimize!(pso, opts)

@test res.fbest <= 1.0e-5

for i in 1:N
    @test res.xbest[i] >= LB[i]
    @test res.xbest[i] <= UB[i]
end

