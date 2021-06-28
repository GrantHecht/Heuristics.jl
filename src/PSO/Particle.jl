
mutable struct Particle{T<:AbstractFloat,N}
    x::SizedArray{Tuple{N},T,1,1,Vector{T}}
    v::SizedArray{Tuple{N},T,1,1,Vector{T}}
    p::SizedArray{Tuple{N},T,1,1,Vector{T}}

    fx::T
    fp::T

    function Particle{T,N}(::UndefInitializer) where {T,N}
        return new{T,N}(SizedVector{N,T}(undef),
                        SizedVector{N,T}(undef),
                        SizedVector{N,T}(undef),
                        T(0.0), T(0.0))
    end
end

# Particle methods
Base.length(p::Particle) = length(p.x)
