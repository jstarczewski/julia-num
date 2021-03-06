using Main.VoronoiDelaunay
using Main.GeometricalPredicates

struct Scaler
    scale::Any
    transx::Any
    transy::Any

    function Scaler(bbox::Array{Int,2})
        _scale = (abs(1 / (bbox[2] - bbox[1])))
        new(_scale, (1 - bbox[1] * _scale), (1 - bbox[1, 2] * _scale))
    end
end

function scaled_point(unscaled_point::Point2D, scaler::Scaler)
    return Point(
        ((getx(unscaled_point) * scaler.scale) + scaler.transx),
        ((gety(unscaled_point) * scaler.scale) + scaler.transy),
    )
end

function unscaled_point(scaled_point::Point2D, scaler::Scaler)
    return Point(
        (getx(scaled_point) - scaler.transx) / scaler.scale,
        (gety(scaled_point) - scaler.transy) / scaler.scale,
    )
end

function scale_p(p, scaler::Scaler)
    return [(p[1] * scaler.scale) + scaler.transx, (p[2] * scaler.scale) + scaler.transy]
end

function unscale_p(p, scaler::Scaler)
    return [(p[1] - scaler.transx) / scaler.scale, (p[2] - scaler.transy) / scaler.scale]
end

function unscale_y(p, scaler)
    return (p - scaler.transx) / scaler.scale
end

function unscale_x(p, scaler)
    return (p - scaler.transy) / scaler.scale
end

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where {T}
    m, n = length(vy), length(vx)
    gx = reshape(repeat(vx, inner = m, outer = 1), m, n)
    gy = reshape(repeat(vy, inner = 1, outer = n), m, n)
    return gx, gy
end

function extractbox(bbox, h0::Real, h1::Real)
    v1 = bbox[1, 1]:h0:bbox[2, 1]
    v2 = bbox[1, 2]:h1:bbox[2, 2]
    return v1, v2
end

function calculateh1(h0::Real)
    return h0 * (sqrt(3) / 2)
end

function shiftevenrows!(x, h0::Real)
    x[2:2:end, :] = x[2:2:end, :] .+ h0 / 2
end

function vectorizededges(triangle)
    a = geta(triangle)
    b = getb(triangle)
    c = getc(triangle)
    ab = Line(Point(getx(a), gety(a)), Point(getx(b), gety(b)))
    ba = Line(Point(getx(b), gety(b)), Point(getx(a), gety(a)))
    bc = Line(Point(getx(b), gety(b)), Point(getx(c), gety(c)))
    cb = Line(Point(getx(c), gety(c)), Point(getx(b), gety(b)))
    ac = Line(Point(getx(a), gety(a)), Point(getx(c), gety(c)))
    ca = Line(Point(getx(c), gety(c)), Point(getx(a), gety(a)))
    return ab, ba, bc, cb, ac, ca
end

function validedges(tess, scaler, fd, geps)
    i = 1
    inside_edges = Array{Line2D,1}()
    outside_edges = Array{Line2D,1}()
    for triangle in tess
        a = geta(triangle)
        b = getb(triangle)
        c = getc(triangle)
        center = unscaled_point(centroid(Primitive(a, b, c)), scaler)
        i += 1
        edges = vectorizededges(triangle)
        if fd([getx(center), gety(center)]) > -geps
            push!(inside_edges, edges...)
        else
            push!(outside_edges, edges...)
        end
    end
    edges = Array{Line2D,1}()
    for edge in delaunayedges(tess)
        push!(edges, Line(geta(edge), getb(edge)))
    end
    inside_edges = filter(edge -> !(edge in outside_edges), inside_edges)
    edges = filter(edge -> !(edge in inside_edges), edges)
    return [edge for edge in edges]
end

function plotedges(edges)
    x = Array{Float64,1}()
    y = Array{Float64,1}()
    for edge in edges
        push!(x, getx(geta(edge)))
        push!(x, getx(getb(edge)))
        push!(x, NaN)
        push!(y, gety(geta(edge)))
        push!(y, gety(getb(edge)))
        push!(y, NaN)
    end
    return x, y
end

function pointstoforces(edges, scaler, Fscale, pfix)
    bars = Array{Point2D,1}()
    barvec = Array{Array{Float64,1},1}()
    points_to_fvces = Dict{Point2D,Array{Float64,1}}()
    for edge in edges
        b = unscaled_point(getb(edge), scaler)
        a = unscaled_point(geta(edge), scaler)
        push!(
            bars,
            Point(getx(a) + ((getx(b) - getx(a)) / 2), gety(a) + ((gety(b) - gety(a)) / 2)),
        )
        push!(
            barvec,
            [
                getx(b) - getx(a)
                gety(b) - gety(a)
            ],
        )
        push!(points_to_fvces, geta(edge) => [0, 0])
        push!(points_to_fvces, getb(edge) => [0, 0])
    end
    L = [sqrt(sum(v_sum .^ 2)) for v_sum in barvec]
    hbars = [fh(getx(p), gety(p)) for p in bars]
    L0 = hbars * Fscale * sqrt(sum(L .^ 2) / sum(hbars .^ 2))
    Fvec = maximum.(L0 - L) ./ L .* barvec
    iterator = 1
    for edge in edges
        a = geta(edge)
        b = getb(edge)
        prev_a = points_to_fvces[a]
        prev_b = points_to_fvces[b]
        push!(points_to_fvces, a => prev_a + (-Fvec[iterator]))
        push!(points_to_fvces, b => prev_b + (Fvec[iterator]))
        iterator = iterator + 1
    end
    for p in pfix
        delete!(points_to_fvces, p)
        push!(points_to_fvces, p => [0, 0])
    end
    return points_to_fvces
end

function finalpositions(points_to_fvces, scaler, deltat, fd, geps, deps, h0)
    new_p = Array{Point2D,1}()
    d_points = Array{Array{Float64,1},1}()
    for (point, force) in points_to_fvces
        p = unscaled_point(point, scaler)
        np = [getx(p), gety(p)] + deltat * force
        push!(new_p, Point2D(np[1], np[2]))
        d = fd(np)
        if d < -geps
            push!(d_points, force)
        end
    end
    final_p = Array{Array{Float64,1},1}()
    for p in new_p
        x = getx(p)
        y = gety(p)
        d = fd([x, y])
        if d > 0
            dgradx = (fd([x + deps, y]) - d) / deps
            dgrady = (fd([x, y + deps]) - d) / deps
            res = [x, y] - [d * dgradx, d * dgrady]
            push!(final_p, res)
        else
            push!(final_p, [x, y])
        end
    end
    final_p = map(p -> scaled_point(Point(p[1], p[2]), scaler), unique(final_p))
    return final_p, moveindex(d_points, deltat, h0)
end

function moveindex(d_points, deltat, h0)
    d = map(row -> sum(deltat * row .^ 2), d_points)
    push!(d, 0)
    return maximum(sqrt.(d) / h0)
end

function movein(p)
    newx = getx(p)
    newy = gety(p)
    if getx(p) <= 1
        newx = round(newx, digits = 8) + 2 *eps(Float64)
    end
    if gety(p) <= 1
        newy = round(newy, digits = 8) + 2 *eps(Float64)
    end
    if getx(p) >= 2
        newx = round(newx, digits = 8) - 3 * eps(Float64)
    end
    if gety(p) >= 2
        newy = round(newy, digits = 8) - 3 * eps(Float64)
    end
    return Point(newx, newy)
end

function distmesh(
    fd,
    fh,
    bbox,
    h0,
    pfix = [],
    ttol = 0.1,
    geps = 0.001 * h0,
    Fscale = 1.2,
    dptol = 0.001,
    deltat = 0.2,
    deps = sqrt(eps(Float64)) * h0,
)
    h1 = calculateh1(h0)
    pold = Inf
    scaler = Scaler(bbox)
    v1, v2 = extractbox(bbox, h0, h1)
    x, y = meshgrid(v1, v2)
    shiftevenrows!(x, h0)
    p = [x[:] y[:]]
    p = [Point(row[1], row[2]) for row in eachrow([x[:] y[:]]) if fd(row) < geps]
    r0 = [1 ./ fh(getx(point), gety(point)) .^ 2 for point in p]
    r0_max = maximum(r0 ./ maximum(r0))
    scaler = Scaler(bbox)
    pfix = [scaled_point(Point(row[1], row[2]), scaler) for row in eachrow(pfix)]
    pfix = map(p -> movein(p), pfix)
    p = [scaled_point(point, scaler) for point in p if (rand(Float64, size(p)))[1] < r0_max]
    p = [p; pfix]
    p = map(p -> movein(p), p)
    p = unique(p)
    while true
        tess = DelaunayTessellation(convex = false)
        push!(tess, p)
        edges = validedges(tess, scaler, fd, geps)
        points_to_fvces = pointstoforces(edges, scaler, Fscale, pfix)
        final_p, move_index =
            finalpositions(points_to_fvces, scaler, deltat, fd, geps, deps, h0)
        if move_index < dptol
            x, y = plotedges(edges)
            break
        else
            p = map(p -> movein(p), final_p)
        end
    end
    x = map(x -> unscale_x(x, scaler), x)
    y = map(y -> unscale_y(y, scaler), y)
    return x, y
end
