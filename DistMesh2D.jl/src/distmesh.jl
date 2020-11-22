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
