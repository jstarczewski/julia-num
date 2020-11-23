using Gadfly
using Main.DistMesh2D

function fh(x, y)::Real
    return 1.0
end

function dc(p)
    return -1 + abs(4 - sqrt(sum(p .^ 2)))
end

function dca(p)
    return sqrt(sum(p .^ 2)) - 2
end

function dcrcc(p)
    return max(drectangle(p, -100.0, 100.0, -100.0, 100.0), - dcircle(p, 20.0, 10.0, 15.0))
end

function eplot(x,y)
    Gadfly.plot(x = x, y = y, Geom.path, Coord.cartesian(fixed = true))
end

#x, y = distmesh2d(dcrcc, fh, [-100.0 -100.0; 100.0 100.0], 7.0, [-100.0 -100.0; -100.0 100.0; 100.0 -100.0; 100.0 100.0])
function fdd(p)
    return dcircle(p, 0.0, 0.0, 1.0)
end
x, y = distmesh2d(fdd, fh, [-1.0 -1.0; 1.0 1.0], 0.2)
eplot(x,y)
