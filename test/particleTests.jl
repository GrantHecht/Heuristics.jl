
using Heuristics, StaticArrays

# Particle State
N = 10
T = typeof(10.0)
p = Heuristics.ParticleState{T,N}(undef)

@test p.x isa SizedArray{Tuple{N},T,1,1}
@test p.v isa SizedArray{Tuple{N},T,1,1}
@test p.p isa SizedArray{Tuple{N},T,1,1}

tempX = rand(N)
tempV = rand(N)
tempP = rand(N)

p.x .= tempX
p.v .= tempV
p.p .= tempP

@test p.x == tempX
@test p.v == tempV
@test p.p == tempP

# Particle Particle Fitness
U = typeof(0.0)
u = Heuristics.ParticleFitness{U}(0,0)

@test u.fx == 0.0
@test u.fp == 0.0
@test u.fx isa U
@test u.fp isa U

u.fx = 5;
u.fp = 6;

@test u.fx == 5.0
@test u.fp == 6.0
@test u.fx isa U
@test u.fp isa U

# Particle Type 
pp = Heuristics.Particle{T,U,N}(undef)

@test pp.state.x isa SizedArray{Tuple{N},T,1,1}
@test pp.state.v isa SizedArray{Tuple{N},T,1,1}
@test pp.state.p isa SizedArray{Tuple{N},T,1,1}
@test pp.fit.fx == 0.0
@test pp.fit.fp == 0.0
@test pp.fit.fx isa U
@test pp.fit.fp isa U

tempX = rand(N)
tempV = rand(N)
tempP = rand(N)

pp.state.x .= tempX
pp.state.v .= tempV
pp.state.p .= tempP
pp.fit.fx = 5;
pp.fit.fp = 6;

@test pp.fit.fx == 5.0
@test pp.fit.fp == 6.0
@test pp.fit.fx isa U
@test pp.fit.fp isa U
@test pp.state.x == tempX
@test pp.state.v == tempV
@test pp.state.p == tempP
@test pp.state.x isa SizedArray{Tuple{N},T,1,1}
@test pp.state.v isa SizedArray{Tuple{N},T,1,1}
@test pp.state.p isa SizedArray{Tuple{N},T,1,1}

