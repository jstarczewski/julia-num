using Test
using LinearAlgebra
using VoronoiDelaunay
using Gadfly
using GeometricalPredicates

# Arguments
# fd - funkcja odległości

struct Scaler
    scale::Any
    transx::Any
    transy::Any

    function Scaler(bbox::Array{Int,2})
        _scale = (abs(1 / (bbox[2] - bbox[1])))
        new(
            round(_scale, digits=4),
            round((1 + (0 - bbox[1] * _scale)), digits=4),
            round((1 + (0 - bbox[1, 2] * _scale)), digits=4),
        )
    end
end


function scaled_point(unscaled_point::Point2D, scaler::Scaler)
    return Point(
        round(((getx(unscaled_point) * scaler.scale) + scaler.transx), digits=4),
        round(((gety(unscaled_point) * scaler.scale) + scaler.transy), digits=4)
    )
end

function unscaled_point(scaled_point::Point2D, scaler::Scaler)
    return Point(
        round((getx(scaled_point) - scaler.transx) / scaler.scale, digits=4),
        round((gety(scaled_point) - scaler.transy) / scaler.scale, digits=4)
    )
end


function fh(x, y)::Real
    #  @info "Value is" s = 0.2
    return 1
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
    pb = [
        Point(row[1], row[2])
        for row in eachrow([x[:] y[:]]) if fd(row) < geps
    ]
    r0 = [1 ./ fh(getx(point), gety(point)) .^ 2 for point in pb]
    r0_max = maximum(r0 ./ maximum(r0))
    pold = Inf
    scaler = Scaler(bbox)
    pps = pb
    pb = [
        scaled_point(point, scaler)
        for point in pb if (rand(Float64, size(pb)))[1] < r0_max
    ]
    while true
        @info "Loop"
        pb = unique!(pb)
        tesselation = DelaunayTessellation(size(pb)[1])
        @info "Push to first tess"
        for p in pb
            @info "P is" p
            push!(tesselation, p)
        end
        interior_points = Array{Point2D,1}()
        for triangle in tesselation
            center = centroid(Primitive(geta(triangle), getb(triangle), getc(triangle)))
            unscaled_p = unscaled_point(Point(getx(center), gety(center)), scaler)
            d = fd([getx(unscaled_p) gety(unscaled_p)])
            if d < -geps
                push!(interior_points, geta(triangle))
                push!(interior_points, getb(triangle))
                push!(interior_points, getc(triangle))
            end
        end
        interior_points = unique!(interior_points)
        interior_tesselation = DelaunayTessellation(size(interior_points)[1])
        @info "Push to second tess without triangles"
        push!(interior_tesselation, interior_points)
        bars = Array{Point2D,1}()
        barvec = Array{Array{Float64,1},1}()
        points_to_fvces = Dict{Point2D,Array{Float64,1}}()
        @info "First iteration over edges started"

        #Point2D(1.4, 1.3464101615137753)
        #Point2D(1.4002786971987105, 1.3468928792218693)
        for edge in delaunayedges(interior_tesselation)
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
        @info "First iteration over edges finished"
        L = [sqrt(sum(v_sum .^ 2)) for v_sum in barvec]
        iterator = 1
        hbars = [fh(getx(p), gety(p)) for p in bars]
        L0 = hbars * Fscale * sqrt(sum(L .^ 2) / sum(hbars .^ 2))
        F = maximum(L0-L)
        @info "Force is " F
        #Fvec=F./L*[1,1].*barvec
        Fvec = F ./ L .* barvec
        @info "Final iteration started"
        for edge in delaunayedges(interior_tesselation)
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
                dgradx = (fd([x + deps, y])-d)
                dgrady = (fd([x, y + deps])-d)
                rdgradx = dgradx/deps
                rdgrady = dgrady/deps
              res = [x, y] - [d*rdgradx, d*rdgrady]
                push!(final_p, res)
            else
                push!(final_p, [x, y])
            end
        end
        pb = map(p -> scaled_point(Point2D(p[1], p[2]), scaler), final_p)
        d_points = map(row -> sum(deltat * row .^ 2), d_points)
        @info "Checking break condition"
        if maximum(sqrt.(d_points) / h0) < dptol
            @info "Breaking out of loop"
            break
        end
   end
   tex = DelaunayTessellation(size(pb)[1])
   @info "Push to output tess"
   @info "Pb size is" size(pb)
   pb = unique!(pb)
   push!(tex, pb)
   x, y = getplotxy(delaunayedges(tex))
   x = map(p -> (p - scaler.transx) / scaler.scale, x)
   y = map(p -> (p - scaler.transy) / scaler.scale, y)
   @info "Plotting started"
   Gadfly.plot(x = x, y=y, Geom.path, Coord.cartesian(fixed=true))
end

function dc(p)
    #return -0.3 + abs(0.7 - sqrt(sum(p .^ 2)))
    return sqrt(sum(p .^ 2)) - 1
end
