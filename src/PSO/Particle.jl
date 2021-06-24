
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

