
using Heuristics
using BenchmarkTools

# Schwefel Function
function schaffer(x)
    obj = 0.5 + (sin(x[1]^2 + x[2]^2)^2 - 0.5)/(1 + 0.001*(x[1]^2+x[2]^2))^2
    return obj 
end

# Setup Problem
d = 2
LB = -100*ones(d)
UB = 100*ones(d)

prob = Problem(schaffer, LB, UB)
pso  = PSO(prob; numParticles = 10)

# optimize
opts = Options(;display = false)
@benchmark optimize!($pso, $opts)