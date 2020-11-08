using Pkg
using Gadfly

try
Pkg.rm("VoronoiDelaunay")
Pkg.rm("GeometricalPredicates")
catch e
    println("Not all packages were installed before")
end
include("GeometricalPredicates.jl/src/GeometricalPredicates.jl")
println("GeometricalPredicates.jl included")

"""
Remember to change GeometricalPredicates module in VoronoiDelaunay.jl/src/VoronoiDelaunay.jl to

using Main.GeometricalPredicates
const GP = Main.GeometricalPredicates
import Main.GeometricalPredicates: geta, getb, getc

"""

include("VoronoiDelaunay.jl/src/VoronoiDelaunay.jl")
println("VoronoiDelaunay.jl included")
println("DistMeshVD.jl included")
include("DistMeshVD.jl")

function fd(p)
     return max(dcb(p, -1,1,-1,1), -dcrc(p, 0,0,0.4))
 end

println("Generating square with hole")
x, y = generate(fd, fh, [-1 -1; 1 1], 0.1, [-1 -1; -1 1; 1 -1; 1 1])
println("Plotting with Gadfly")
Gadfly.plot(x = x, y = y, Geom.path, Coord.cartesian(fixed = true))
println("Generation compleated")
