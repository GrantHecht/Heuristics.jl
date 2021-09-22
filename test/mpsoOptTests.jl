using Heuristics, Suppressor 

include("testProblems.jl")

begin
    N = 14
    M = 200
    LB = -5 .* ones(N)
    UB = 5 .* ones(N)
    prob1 = Problem(sphereFunc, LB, UB)
    opts = Options(;display = false, maxStallIters = 1000)
    mpso1 = MPSO(prob1; numParticles = M)

    res = optimize!(mpso1, opts)

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
    opts = Options(;display = true, maxStallIters = 1000)
    mpso2 = MPSO(prob2; numParticles = M, MFDType = :Average)

    @suppress_out begin
        optimize!(mpso2, opts)
    end

    @test mpso2.results.fbest <= 1.0e-5

    for i in 1:N
        @test mpso2.results.xbest[i] >= LB[i]
        @test mpso2.results.xbest[i] <= UB[i]
    end
end
 
# Sphere function test 3
begin
    N = 2
    M = 1000
    prob3 = Problem(sphereFunc, N)
    opts = Options(;display = true, maxStallIters = 1000)
    mpso3 = MPSO(prob3; numParticles = M, MFDType = :Global)

    @suppress_out begin
        optimize!(mpso3, opts)
    end

    @test mpso3.results.fbest <= 1.0e-5

    for i in 1:N
        @test mpso3.results.xbest[i] >= LB[i]
        @test mpso3.results.xbest[i] <= UB[i]
    end
end
 