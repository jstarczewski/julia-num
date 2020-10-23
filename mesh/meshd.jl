using Deldir
using Gadfly
using GeometricalPredicates

struct Scaler
    scale::Any
    transx::Any
    transy::Any

    function Scaler(bbox::Array{Int,2})
        _scale = (abs(1 / (bbox[2] - bbox[1])))
        new(_scale, (0 - bbox[1] * _scale), (0 - bbox[1, 2] * _scale))
    end
end


function scale_p(p, scaler::Scaler)
    return [(p[1] * scaler.scale) + scaler.transx, (p[2] * scaler.scale) + scaler.transy]
end

function unscale_p(p, scaler::Scaler)
    return [(p[1] - scaler.transx) / scaler.scale, (p[2] - scaler.transy) / scaler.scale]
end

function scaled_point(unscaled_point::Point2D, scaler::Scaler)
    return Point(
        round(((getx(unscaled_point) * scaler.scale) + scaler.transx), digits = 4),
        round(((gety(unscaled_point) * scaler.scale) + scaler.transy), digits = 4),
    )
end

function unscaled_point(scaled_point::Point2D, scaler::Scaler)
    return Point(
        round((getx(scaled_point) - scaler.transx) / scaler.scale, digits = 4),
        round((gety(scaled_point) - scaler.transy) / scaler.scale, digits = 4),
    )
end

function dc(p)
    #return -0.3 + abs(0.7 - sqrt(sum(p .^ 2)))
    return sqrt(sum(p .^ 2)) - 1
end

function fh(x, y)::Real
    #  @info "Value is" s = 0.2
    return 0.2
end

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where {T}
    m, n = length(vy), length(vx)
    #  @info "Length m = " m
    #  @info "Length n = " n
    gx = reshape(repeat(vx, inner = m, outer = 1), m, n)
    gy = reshape(repeat(vy, inner = 1, outer = n), m, n)
    #  @info "Grid x" gx
    #  @info "Grid y" gy
    return gx, gy
end

#bbox = [0 0; 1 1]

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
    final_p = Array{Array{Float64,1},1}()
    shiftevenrows!(x, h0)
    po = [x[:] y[:]]
    po = [vcat(row) for row in eachrow(po) if fd(row) < -geps]
    r0 = [1 ./ fh(row...) .^ 2 for row in po]
    r0_max = maximum(r0 ./ maximum(r0))
    pold = Inf
    scaler = Scaler(bbox)
    po = [vcat(scale_p(row, scaler)) for row in po if (rand(Float64, size(po)))[1] < r0_max]
    pa = transpose(reshape(vcat(po...), 2, length(po)))
    del, vor, summ = deldir(pa[:, 1], pa[:, 2])
    #while true
        println(typeof(pa[:,1]))
        println(typeof(pa[:,2]))
        Dx, Dy = edges(del)
        interior_points = Array{Point2D,1}()
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
        for tp in tri_point
            center = unscaled_point(centroid(Primitive(tp[1], tp[2], tp[3])), scaler)
            if fd([getx(center), gety(center)]) < -geps
                push!(interior_points, unscaled_point(tp[1], scaler))
                push!(interior_points, unscaled_point(tp[2], scaler))
                push!(interior_points, unscaled_point(tp[3], scaler))
            end
        end
        interior_points =
            map(p -> scale_p([getx(p), gety(p)], scaler), unique!(interior_points))
        interior_points =
            transpose(reshape(vcat(interior_points...), 2, length(interior_points)))
        del, vor, summ = deldir(interior_points[:, 1], interior_points[:, 2])
        dedges = [Line(unscaled_point(Point(r[1], r[2]), scaler), unscaled_point(Point(r[3], r[4]), scaler)) for r in eachrow(del)]
        bars = Array{Point2D,1}()
        barvec = Array{Array{Float64,1},1}()
        points_to_fvces = Dict{Point2D,Array{Float64,1}}()
        for edge in dedges
            b = unscaled_point(getb(edge), scaler)
            a = unscaled_point(geta(edge), scaler)
            # Obliczamy jako wartość fh() w połowie kazdego bara
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
            push!(points_to_fvces, geta(edge) => [0, 0])
            push!(points_to_fvces, getb(edge) => [0, 0])
        end
        L = [sqrt(sum(v_sum .^ 2)) for v_sum in barvec]
        iterator = 1
        hbars = [fh(getx(p), gety(p)) for p in bars]
        L0 = hbars * Fscale * sqrt(sum(L .^ 2) / sum(hbars .^ 2))
        sqrt(sum(L .^ 2) / sum(hbars .^ 2))
        F = maximum(L0 - L)
        @info "Force is " F
        #Fvec=F./L*[1,1].*barvec
        Fvec = F ./ L .* barvec
        @info "Final iteration started"
        for edge in dedges
            prev_a = points_to_fvces[geta(edge)]
            prev_b = points_to_fvces[getb(edge)]
            push!(points_to_fvces, geta(edge) => prev_a + (-Fvec[iterator]))
            push!(points_to_fvces, getb(edge) => prev_b + (Fvec[iterator]))
            iterator = iterator + 1
        end
        @info "Final iteration finished"
        new_p = Array{Point2D,1}()
        d_points = Array{Array{Float64,1},1}()
        for (point, force) in points_to_fvces
            p = unscaled_point(point, scaler)
            np = [getx(p), gety(p)] + 0.2 * force
            push!(new_p, Point2D(np[1], np[2]))
            d = fd(np)
            if d < -geps
                push!(d_points, force)
            end
        end
        @info "Mapping points"
        new_p_raw = map(p -> [getx(p), gety(p)], new_p)
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
        final_p = map(p -> scaled_point(Point2D(p[1], p[2]), scaler), final_p)
        println(final_p)
        final_p = map(p -> [getx(p), gety(p)], final_p)
        #pa = [vcat(scale_p([round(p[1],RoundNearestTiesAway, sigdigits=4), round(p[2], RoundNearestTiesAway, sigdigits=4)], scaler)) for p in final_p]
        pa = transpose(reshape(vcat(final_p...), 2, length(final_p)))
        println(pa)
        d_points = map(row -> sum(deltat * row .^ 2), d_points)
        @info "Checking break condition"
        if maximum(sqrt.(d_points) / h0) < dptol
            @info "Breaking out of loop"
            #break
        end
    #end
    del, vor, summ = deldir(pa[:,1], pa[:,2])

    Dx, Dy = edges(del)
    Gadfly.plot(x=Dx, y=Dy, Geom.path, Coord.cartesian(fixed=true))
end
