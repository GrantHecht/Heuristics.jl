
struct ParticleState{T<:AbstractFloat, N}
    x::SizedArray{Tuple{N},T,1,1,Vector{T}}
    v::SizedArray{Tuple{N},T,1,1,Vector{T}}
    p::SizedArray{Tuple{N},T,1,1,Vector{T}}

    function ParticleState{T,N}(::UndefInitializer) where {T,N}
        return new{T,N}(Vector{T}(undef, N), 
                        Vector{T}(undef, N), 
                        Vector{T}(undef, N))
    end
end

mutable struct ParticleFitness{U<:Real}
    fx::U
    fp::U 

    function ParticleFitness{U}(fxVal::Real,fpVal::Real) where {U<:Real}
        return new{U}(U(fxVal), U(fpVal))
    end
end

struct Particle{T<:AbstractFloat, U<:Real, N}
    state::ParticleState{T,N}
    fit::ParticleFitness{U}

    function Particle{T,U,N}(::UndefInitializer) where {T,U,N}
        return new{T,U,N}(ParticleState{T,N}(undef),
                          ParticleFitness{U}(0.0, 0.0))
    end
end

