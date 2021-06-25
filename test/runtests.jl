using Heuristics
using Test, SafeTestsets

@time begin
@time @safetestset "Particle tests..." begin include("particleTests.jl") end     
@time @safetestset "Swarm tests..." begin include("swarmTests.jl") end
end