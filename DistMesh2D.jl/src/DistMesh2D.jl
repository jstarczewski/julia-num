module DistMesh2D

using GeometricalPredicates, Deldir

include("scaler.jl")
include("util.jl")
include("distancefunctions.jl")
include("distmesh.jl")

export distmesh2d, drectangle, dcircle

end
