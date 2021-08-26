using Heuristics
using Test, SafeTestsets

@time begin
#@time @safetestset "Particle tests..." begin include("particleTests.jl") end     
#@time @safetestset "Swarm tests..." begin include("swarmTests.jl") end
@time @safetestset "PSO optimization test..." begin include("psoOptTests.jl") end
end