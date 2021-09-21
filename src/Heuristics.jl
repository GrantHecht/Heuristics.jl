module Heuristics

using StaticArrays
using Format 
using ThreadPools
import Random: shuffle!

# Required
include("Problem.jl")
include("Options.jl")
include("Optimizers.jl")
include("Results.jl")

# PSO includes
include("./PSO/psoUtil.jl")
include("./PSO/Particle.jl")
include("./PSO/Swarm.jl")
include("./PSO/PSO.jl")
include("./PSO/MS_PSO.jl")
include("./PSO/MPSO.jl")

# Exports
export Problem
export Options
export optimize!

export PSO
export MS_PSO
export MPSO

end
