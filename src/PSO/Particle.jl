
mutable struct Particle{T<:AbstractFloat, U<:Real, N}
    x::SizedArray{Tuple{N},T,1,1,Vector{T}}
    v::SizedArray{Tuple{N},T,1,1,Vector{T}}
    p::SizedArray{Tuple{N},T,1,1,Vector{T}}

    fx::U
    fp::U

    function Particle{T,U,N}(::UndefInitializer) where {T,U,N}
        return new{T,U,N}(SizedVector{N,T}(undef),
                          SizedVector{N,T}(undef),
                          SizedVector{N,T}(undef),
                          U(0.0), U(0.0))
    end
end

# Particle methods
Base.length(p::Particle) = length(p.x)
