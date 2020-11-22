function distmesh2d(
    fd::Function,
    fh::Function,
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
    p = [vcat(row) for row in eachrow(p) if fd(row) < -geps]
    r0 = [1 ./ fh(row...) .^ 2 for row in p]
    r0max = maximum(r0 ./ maximum(r0))
    p = [vcat(scaledpoint(row, scaler)) for row in p if (rand(Float64, size(p)))[1] < r0max]
    p = transpose(reshape(vcat(p...), 2, length(p)))
    pfix = [vcat(scaledpoint(row, scaler)) for row in eachrow(pfix)]
    pfix = transpose(reshape(vcat(pfix...), 2, length(pfix)))
    while true
        del, vor, summ = deldir(p[:, 1], p[:, 2])
        trigs = triangles(del, summ)
        line_edges = validedges(trigs, del, scaler, fd, geps)
        pointstofvces = pointstoforces(line_edges, scaler, Fscale, pfix, fh)
        finalp, moveindex =
            finalpositions(pointstofvces, scaler, deltat, fd, geps, deps, h0)
        if moveindex < dptol
            x, y = plotedges(line_edges)
            break
        else
            p = finalp
        end
    end
    x = map(x -> unscalex(x, scaler), x)
    y = map(y -> unscaley(y, scaler), y)
    return x, y
end

function calculateh1(h0::Real)
    return h0 * (sqrt(3) / 2)
end

function extractbox(bbox, h0::Real, h1::Real)
    v1 = bbox[1, 1]:h0:bbox[2, 1]
    v2 = bbox[1, 2]:h1:bbox[2, 2]
    return v1, v2
end

function shiftevenrows!(x, h0::Real)
    x[2:2:end, :] = x[2:2:end, :] .+ h0 / 2
end

function vectorizededges(triangle)
    a = triangle[1]
    b = triangle[2]
    c = triangle[3]
    ab = [getx(a) gety(a); getx(b) gety(b)]
    ba = [getx(b) gety(b); getx(a) gety(a)]
    bc = [getx(b) gety(b); getx(c) gety(c)]
    cb = [getx(c) gety(c); getx(b) gety(b)]
    ac = [getx(a) gety(a); getx(c) gety(c)]
    ca = [getx(c) gety(c); getx(a) gety(a)]
    return ab, ba, bc, cb, ac, ca
end

function validedges(triangles, del, scaler, fd, geps)
    i = 1
    inside_edges = Array{Array{Float64,2},1}()
    outside_edges = Array{Array{Float64,2},1}()
    for triangle in triangles
        center = unscaledpoint2d(
            centroid(Primitive(triangle[1], triangle[2], triangle[3])),
            scaler,
        )
        i += 1
        a = triangle[1]
        b = triangle[2]
        c = triangle[3]
        edges = vectorizededges(triangle)
        if fd([getx(center), gety(center)]) > -geps
            push!(inside_edges, edges...)
        else
            push!(outside_edges, edges...)
        end
    end
    nedges = inside_edges
    yedges = outside_edges
    inside_edges = filter(edge -> !(edge in outside_edges), inside_edges)
    edges = [[r[1] r[2]; r[3] r[4]] for r in eachrow(del)]
    edges = filter(edge -> !(edge in inside_edges), edges)
    x = Array{Float64,1}()
    y = Array{Float64,1}()
    return [Line(Point(edge[1], edge[3]), Point(edge[2], edge[4])) for edge in edges]
end

function pointstoforces(edges, scaler, Fscale, pfix, fh)
    bars = Array{Point2D,1}()
    barvec = Array{Array{Float64,1},1}()
    points_to_fvces = Dict{Point2D,Array{Float64,1}}()
    for edge in edges
        b = unscaledpoint2d(getb(edge), scaler)
        a = unscaledpoint2d(geta(edge), scaler)
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
    for p in eachrow(pfix)
        p = Point(p[1], p[2])
        delete!(points_to_fvces, p)
        push!(points_to_fvces, p => [0, 0])
    end
    return points_to_fvces
end

function finalpositions(points_to_fvces, scaler, deltat, fd, geps, deps, h0)
    new_p = Array{Point2D,1}()
    d_points = Array{Array{Float64,1},1}()
    for (point, force) in points_to_fvces
        p = unscaledpoint2d(point, scaler)
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
    final_p = unique(final_p)
    final_p = map(p -> scaledpoint([p[1], p[2]], scaler), final_p)
    final_p = transpose(reshape(vcat(final_p...), 2, length(final_p)))
    return final_p, move_index(d_points, deltat, h0)
end

function move_index(d_points, deltat, h0)
    d = map(row -> sum(deltat * row .^ 2), d_points)
    push!(d, 0)
    return maximum(sqrt.(d) / h0)
end
