using Heuristics, StaticArrays

# Particle State
N = 5
M = 5
T = typeof(10.0)
U = typeof(10.0)
V = typeof(0)
s = Heuristics.Swarm{T,U,V,N,M}(undef)

# Initialization
@test s.particles isa SizedArray{Tuple{M},Heuristics.Particle{T,U,N},1,1,Vector{Heuristics.Particle{T,U,N}}}

# Initializing temp variables 
tempXs  = rand(N,M)
tempPs  = rand(N,M)
tempVs  = rand(N,M)
tempD   = rand(N)
tempFxs = rand(M)
tempFps = rand(M)
tempN   = 10.0
tempB   = rand()
tempW   = rand()
tempC   = 11.0
tempY₁  = rand()
tempY₂  = rand()

# Setting values
s.b  = tempB
s.N  = tempN
s.W  = tempW
s.c  = tempC
s.y₁ = tempY₁
s.y₂ = tempY₂
s.d .= tempD
@inbounds for i in 1:M
    s[i].fx = tempFxs[i]
    s[i].fp = tempFps[i]
    for j in 1:N
        s[i].x[j] = tempXs[j]
        s[i].p[j] = tempPs[j]
        s[i].v[j] = tempVs[j]
    end
end

# Testing values
@test s.particles isa SizedArray{Tuple{M},Heuristics.Particle{T,U,N},1,1,Vector{Heuristics.Particle{T,U,N}}}
@test s.b isa U
@test s.d isa SizedArray{Tuple{N},T,1,1,Vector{T}}
@test s.N isa V
@test s.W isa T
@test s.c isa V
@test s.y₁ isa T 
@test s.y₂ isa T 
@test s.b  == U(tempB)
@test s.N  == V(tempN)
@test s.W  == T(tempW)
@test s.c  == V(tempC)
@test s.y₁ == T(tempY₁)
@test s.y₂ == T(tempY₂)
@inbounds for i in 1:N
    @test s.d[i] == T(tempD[i])
end
@inbounds for i in 1:M
    @test s[i].fx == U(tempFxs[i])
    @test s[i].fp == U(tempFps[i])
    for j in 1:N
        @test s[i].x[j] == T(tempXs[j])
        @test s[i].p[j] == T(tempPs[j])
        @test s[i].v[j] == T(tempVs[j])
    end
end
