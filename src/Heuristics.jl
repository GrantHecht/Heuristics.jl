module Heuristics

using StaticArrays
using Format 
import Random: shuffle!

# Required
include("Problem.jl")
include("Options.jl")
include("Optimizers.jl")
include("Results.jl")

# PSO includes
include("./PSO/Particle.jl")
include("./PSO/Swarm.jl")
include("./PSO/PSO.jl")

# Exports
export Problem
export Options
export optimize!

export PSO

end
