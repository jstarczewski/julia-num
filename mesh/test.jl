using Gadfly

function fh(x, y)::Real
    return 1
end

function dc(p)
    return -1 + abs(4 - sqrt(sum(p .^ 2)))
end

function dca(p)
    return sqrt(sum(p .^ 2)) - 2
end

function dcb(p, x1, x2, y1, y2)
    return -min(min(min(-y1 + p[2], y2 - p[2], -x1 + p[1], x2 - p[1])))
end

function dcrc(p, xc, yc, r)
    return sqrt((p[1] - xc) .^ 2 + (p[2] - yc) .^ 2) - r
end

function dcrcc(p)
    return max(dcb(p, -100, 100, -100, 100), -dcrc(p, 20, 10, 15))
end

function eplot(x,y)
    Gadfly.plot(x = x, y = y, Geom.path, Coord.cartesian(fixed = true))
end

x, y = generate(dcrcc, fh, [-100 -100; 100 100], 7.0, [-100 -100; -100 100; 100 -100; 100 100])
