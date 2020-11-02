using Deldir
using Gadfly
using GeometricalPredicates
using Main.PointScaler

function fh(x, y)::Real
    return 1
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

function buildTriangles(del, summ)
    generators = Dict{Int64,Array{Int,1}}()
        for (index, row) in enumerate(eachrow(summ))
            generators[index] = []
        end
        for row in eachrow(del)
            generators[row[5]] = append!(generators[row[5]], row[6])
            generators[row[6]] = append!(generators[row[6]], row[5])
        end
        triangles = Array{Array{Int,1},1}()
        for k in keys(generators)
            for e in generators[k]
                i = findall(in(generators[e]), generators[k])
                common = generators[k][i]
                for c in common
                    tri = sort([c, k, e])
                    if !(tri in triangles)
                        push!(triangles, tri)
                    end
                end
            end
        end
        tri_point = map(
            t -> [
                Point(summ[t[1], 1], summ[t[1], 2]),
                Point(summ[t[2], 1], summ[t[2], 2]),
                Point(summ[t[3], 1], summ[t[3], 2]),
            ],
            triangles,
        )
        return tri_point
end

function buildEdgesToPlot(triangles, del, scaler, fd, geps)
        i = 1
        nedges = Array{Array{Float64,2},1}()
        yedges = Array{Array{Float64,2},1}()
        for tp in triangles
            center = unscaled_point(centroid(Primitive(tp[1], tp[2], tp[3])), scaler)
            i += 1
            if fd([getx(center), gety(center)]) > -geps
                a = tp[1]
                b = tp[2]
                c = tp[3]
                ab = [getx(a) gety(a); getx(b) gety(b)]
                ba = [getx(b) gety(b); getx(a) gety(a)]
                bc = [getx(b) gety(b); getx(c) gety(c)]
                cb = [getx(c) gety(c); getx(b) gety(b)]
                ac = [getx(a) gety(a); getx(c) gety(c)]
                ca = [getx(c) gety(c); getx(a) gety(a)]
                push!(nedges, ab)
                push!(nedges, ba)
                push!(nedges, bc)
                push!(nedges, cb)
                push!(nedges, ac)
                push!(nedges, ca)
            else
                a = tp[1]
                b = tp[2]
                c = tp[3]
                ab = [getx(a) gety(a); getx(b) gety(b)]
                ba = [getx(b) gety(b); getx(a) gety(a)]
                bc = [getx(b) gety(b); getx(c) gety(c)]
                cb = [getx(c) gety(c); getx(b) gety(b)]
                ac = [getx(a) gety(a); getx(c) gety(c)]
                ca = [getx(c) gety(c); getx(a) gety(a)]
                push!(yedges, ab)
                push!(yedges, ba)
                push!(yedges, bc)
                push!(yedges, cb)
                push!(yedges, ac)
                push!(yedges, ca)
            end
        end
        nedges = filter(a -> !(a in yedges), nedges)
        tredges = [[r[1] r[2]; r[3] r[4]] for r in eachrow(del)]
        tredges = filter(a -> !(a in nedges), tredges)
        ppss = Array{Float64,1}()
        ppyy = Array{Float64,1}()
        sedges = Array{Line2D{Point2D},1}()
        for tredge in tredges
            push!(ppss, tredge[1])
            push!(ppss, tredge[2])
            push!(ppss, NaN)
            push!(ppyy, tredge[3])
            push!(ppyy, tredge[4])
            push!(ppyy, NaN)
            push!(sedges, Line(Point(tredge[1], tredge[3]), Point(tredge[2], tredge[4])))
        end
        return sedges, ppss, ppyy
end

function buildForces(sedges, scaler, Fscale)
        bars = Array{Point2D,1}()
        barvec = Array{Array{Float64,1},1}()
        points_to_fvces = Dict{Point2D,Array{Float64,1}}()
        for edge in sedges
            b = unscaled_point(getb(edge), scaler)
            a = unscaled_point(geta(edge), scaler)
            push!(
                bars,
                Point(
                    getx(a) + ((getx(b) - getx(a)) / 2),
                    gety(a) + ((gety(b) - gety(a)) / 2),
                ),
            )
            push!(
                barvec,
                [
                    getx(b) - getx(a)
                    gety(b) - gety(a)
                ],
            )
            if !haskey(points_to_fvces, geta(edge))
                push!(points_to_fvces, geta(edge) => [0, 0])
            end
            if !haskey(points_to_fvces, getb(edge))
                push!(points_to_fvces, getb(edge) => [0, 0])
            end
        end
        L = [sqrt(sum(v_sum .^ 2)) for v_sum in barvec]
        iterator = 1
        hbars = [fh(getx(p), gety(p)) for p in bars]
        L0 = hbars * Fscale * sqrt(sum(L .^ 2) / sum(hbars .^ 2))
        sqrt(sum(L .^ 2) / sum(hbars .^ 2))
        F = maximum.(L0 - L)
        Fvec = F ./ L .* barvec
        for edge in sedges
            prev_a = points_to_fvces[geta(edge)]
            prev_b = points_to_fvces[getb(edge)]
            push!(points_to_fvces, geta(edge) => prev_a + (-Fvec[iterator]))
            push!(points_to_fvces, getb(edge) => prev_b + (Fvec[iterator]))
            iterator = iterator + 1
        end
        return points_to_fvces
end

function applyForcesToPoints(points_to_fvces, scaler, deltat, fd, geps)
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
        return d_points, new_p
end

function buildFinalPoints(d_points, new_p, scaler, deps, fd)
    final_p = Array{Array{Float64,1},1}()
        for p in new_p
            x = getx(p)
            y = gety(p)
            d = fd([x, y])
            if d > 0
                dgradx = (fd([x + deps, y]) - d)
                dgrady = (fd([x, y + deps]) - d)
                rdgradx = dgradx / deps
                rdgrady = dgrady / deps
                res = [x, y] - [d * rdgradx, d * rdgrady]
                push!(final_p, res)
            else
                push!(final_p, [x, y])
            end
    end
    return unique(final_p)
end

function buildMoveIndex(d_points, deltat, h0)
        d = map(row -> sum(deltat * row .^ 2), d_points)
        push!(d, 0)
        return maximum(sqrt.(d) / h0)
end

function generate(fd, fh, bbox, h0)
    ttol = 0.1
    geps = 0.001 * h0
    Fscale = 1.2
    dptol = 0.001
    deps = sqrt(eps(Float64)) * h0
    h1 = calculateh1(h0)
    deltat = 0.2
    v1, v2 = extractbox(bbox, h0, h1)
    x, y = meshgrid(v1, v2)
    shiftevenrows!(x, h0)
    p = [x[:] y[:]]
    p = [vcat(row) for row in eachrow(p) if fd(row) < -geps]
    r0 = [1 ./ fh(row...) .^ 2 for row in p]
    r0_max = maximum(r0 ./ maximum(r0))
    pold = Inf
    scaler = Scaler(bbox)
    p = [vcat(scale_p(row, scaler)) for row in p if (rand(Float64, size(p)))[1] < r0_max]
    p = transpose(reshape(vcat(p...), 2, length(p)))
    ppss = []
    ppyy = []
    while true
        del, vor, summ = deldir(p[:, 1], p[:, 2])
        triangles = buildTriangles(del, summ)
        sedges, ppss, ppyy = buildEdgesToPlot(triangles, del, scaler, fd, geps)
        points_to_fvces = buildForces(sedges, scaler, Fscale)
        d_points, new_p = applyForcesToPoints(points_to_fvces, scaler, deltat, fd, geps)
        moveIndex = buildMoveIndex(d_points, deltat, h0)
        final_p = buildFinalPoints(d_points, new_p, scaler, deps, fd)
        final_p = map(p -> scale_p([p[1], p[2]], scaler), final_p)
        final_p = transpose(reshape(vcat(final_p...), 2, length(final_p)))
        if moveIndex < dptol
            break
        else
            p = final_p
        end
    end
    Gadfly.plot(x = ppss, y = ppyy, Geom.path, Coord.cartesian(fixed = true))
end
