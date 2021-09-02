using Heuristics, Suppressor

include("testProblems.jl")

# Sphere function test 1
begin
    N = 14
    M = 200
    LB = -5 .* ones(N)
    UB = 5 .* ones(N)
    prob1 = Problem(sphereFunc, LB, UB)
    opts = Options(;display = false, useParallel=true)
    pso1 = MS_PSO(prob1; numParticlesPerSwarm = M)

    res = optimize!(pso1, opts)

    @test res.fbest <= 1.0e-5

    for i in 1:N
        @test res.xbest[i] >= LB[i]
        @test res.xbest[i] <= UB[i]
    end
end

# Sphere function test 2
begin
    N = 2
    M = 1000
    prob2 = Problem(sphereFunc, N)
    opts = Options(;display = true, useParallel = true)
    pso2 = MS_PSO(prob2; numParticlesPerSwarm = M, initMethod = :LogisticsMap)

    @suppress_out begin
        optimize!(pso2, opts)
    end

    @test pso2.results.fbest <= 1.0e-5

    for i in 1:N
        @test pso2.results.xbest[i] >= LB[i]
        @test pso2.results.xbest[i] <= UB[i]
    end
end

