using Heuristics, Suppressor

include("testProblems.jl")

# Sphere function test 1
begin
    N = 14
    M = 200
    LB = -5 .* ones(N)
    UB = 5 .* ones(N)
    prob1 = Problem{N}(sphereFunc, LB, UB)
    opts = Options(;display = false, useParallel=true)
    pso1 = PSO{M}(prob1)

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
    prob2 = Problem{N}(sphereFunc)
    pso2 = PSO{M}(prob2; initMethod = :LogisticsMap)

    @suppress_out begin
        optimize!(pso2)
    end

    @test pso2.results.fbest <= 1.0e-5

    for i in 1:N
        @test pso2.results.xbest[i] >= LB[i]
        @test pso2.results.xbest[i] <= UB[i]
    end
end

