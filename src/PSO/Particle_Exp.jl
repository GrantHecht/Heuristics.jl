
# Experimental perticle class based on StaticArray's SizedArray type UNFINISHED

struct ParticleStates{S<:Tuple,T<:AbstractFloat,N,M,TData<:AbstractArray{T,M}} <: StaticArray{S,T,N} 
    x::TData # Current Position
    v::TData # Current Velocity
    p::TData # Best Position

    function ParticleStates{S,T,N,M,TData}(a::TData) where {S,T,N,M,TData<:AbstractArray{T,M}}
        StaticArrays.require_one_based_indexing(a)
        if size(a) != size_to_tuple(S) && size(a) != (tuple_prod(S),)
            throw(DimensionMismatch("Dimensions $(size(a)) don't match static size $S"))
        end
        return new{S,T,N,M,TData}(a,deepcopy(a),deepcopy(a))
    end

    function ParticleStates{S,T,N,1,TData}(::UndefInitializer) where {S,T,N,TData<:AbstractArray{T,1}}
        return new{S,T,N,1,TData}(TData(undef, StaticArrays.tuple_prod(S)), 
                                  TData(undef, StaticArrays.tuple_prod(S)),
                                  TData(undef, StaticArrays.tuple_prod(S)))
    end

    function ParticleStates{S,T,N,N,TData}(::UndefInitializer) where {S,T,N,TData<:AbstractArray{T,N}}
        return new{S,T,N,N,TData}(TData(undef, StaticArrays.size_to_tuple(S)...),
                                  TData(undef, StaticArrays.size_to_tuple(S)...),
                                  TData(undef, StaticArrays.size_to_tuple(S)...))
    end 
end

function ParticleStates{S,T,N,N}(::UndefInitializer) where {S,T,N}
    return ParticleStates{S,T,N,N,Array{T,N}}(undef)
end
function ParticleStates{S,T,N,1}(::UndefInitializer) where {S,T,N}
    return ParticleStates{S,T,N,1,Vector{T}}(undef)
end
@inline function ParticleStates{S,T,N}(::UndefInitializer) where {S,T,N}
    return ParticleStates{S,T,N,N}(undef)
end 
@inline function ParticleStates{S,T}(::UndefInitializer) where {S,T}
    return ParticleStates{S,T,StaticArrays.tuple_length(S)}(undef) 
end



#mutable struct ParticleFitness{}
#    fx::
#    fp::
#end