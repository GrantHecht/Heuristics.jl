
mutable struct Particle{T<:AbstractFloat}
    x::Vector{T}
    v::Vector{T}
    p::Vector{T}

    fx::T
    fp::T

    function Particle{T}(nDims::Integer) where {T}
        if nDims < 0
            throw(ArgumentError("nDims cannot be less than 0."))
        end
        return new{T}(Vector{T}(undef, nDims),
                      Vector{T}(undef, nDims),
                      Vector{T}(undef, nDims),
                      T(0.0), T(0.0))
    end
 end

function Particle{T}(::UndefInitializer) where {T}
    return Particle{T}(Vector{T}(undef, 0),
                    Vector{T}(undef, 0),
                    Vector{T}(undef, 0),
                    T(0.0), T(0.0))
end

# Particle methods
Base.length(p::Particle) = length(p.x)
