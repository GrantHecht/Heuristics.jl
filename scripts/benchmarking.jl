
using Heuristics
using BenchmarkTools

# Schwefel Function
function schaffer(x)
    obj = 0.5 + (sin(x[1]^2 + x[2]^2)^2 - 0.5)/(1 + 0.001*(x[1]^2+x[2]^2))^2
    return obj 
end

function waveDrop(x)
    obj = -(1 + cos(12*sqrt(x[1]^2 + x[2]^2)))/(0.5*(x[1]^2 + x[2]^2) + 2.0)
    return obj
end

# Setup Problem
d = 2
LB = -5.12*ones(d)
UB = 5.12*ones(d)

prob = Problem(waveDrop, LB, UB)
pso  = PSO(prob; numParticles = 100)

# optimize
opts = Options(;display = false, maxStallIters = 100)
@benchmark optimize!($pso, $opts)
#optimize!(pso, opts)