
using Heuristics, StaticArrays

# Particle Type 
N = 10
T = typeof(10.0)
U = typeof(10.0)
pp = Heuristics.Particle{T,U,N}(undef)

@test pp.x isa SizedArray{Tuple{N},T,1,1}
@test pp.v isa SizedArray{Tuple{N},T,1,1}
@test pp.p isa SizedArray{Tuple{N},T,1,1}
@test pp.fx == 0.0
@test pp.fp == 0.0
@test pp.fx isa U
@test pp.fp isa U

tempX = rand(N)
tempV = rand(N)
tempP = rand(N)

pp.x .= tempX
pp.v .= tempV
pp.p .= tempP
pp.fx = 5;
pp.fp = 6;

@test pp.fx == 5.0
@test pp.fp == 6.0
@test pp.fx isa U
@test pp.fp isa U
@test pp.x == tempX
@test pp.v == tempV
@test pp.p == tempP
@test pp.x isa SizedArray{Tuple{N},T,1,1}
@test pp.v isa SizedArray{Tuple{N},T,1,1}
@test pp.p isa SizedArray{Tuple{N},T,1,1}

