mutable struct Swarm{T<:AbstractFloat,U<:Real,V<:Int,N,M}

    particles::SizedArray{Tuple{M},Particle{T,U,N},1,1,Vector{Particle{T,U,N}}}

    # Global best objective function value
    b::U

    # Location of global best objective function value 
    d::SizedArray{Tuple{N},T,1,1,Vector{T}}

    n::V    # Neighborhood size 
    w::T    # Inertia
    c::V    # Adaptive inertia counter

    y₁::T   # Self adjustment weight
    y₂::T   # Social adjustment weight

    function Swarm{T,U,V,N,M}(::UndefInitializer) where {T,U,V,N,M}
        return new{T,U,V,N,M}(
            SizedVector{M}([Particle{T,U,N}(undef) for i in 1:M]),
            U(0.0), SizedVector{N,T}(undef), 
            V(0), T(0.0), V(0), T(0.0), T(0.0))
    end
end

# Swarm methods
Base.length(s::Swarm) = length(s.particles)

function Base.getindex(s::Swarm, i::Int)
    1 <= i <= length(s.particles) || throw(BoundsError(s, i))
    return s.particles[i]
end

function Base.setindex!(s::Swarm, v, i::Int)
    1 <= i <= length(s.particles) || throw(BoundsError(s, i))
    s.particles[i] = v
end


