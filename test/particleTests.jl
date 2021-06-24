
using Heuristics, StaticArrays

N = 10
T = typeof(10.0)
p = Heuristics.ParticleState{T,N}(undef)

@test p.x isa SizedArray{Tuple{N},T,1,1}
@test p.v isa SizedArray{Tuple{N},T,1,1}
@test p.p isa SizedArray{Tuple{N},T,1,1}

tempX = rand(10)
tempV = rand(10)
tempP = rand(10)

p.x .= tempX
p.v .= tempV
p.p .= tempP

@test p.x == tempX
@test p.v == tempV
@test p.p == tempP