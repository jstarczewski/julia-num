using Test
using Gadfly
using LinearAlgebra
using VoronoiDelaunay
using GeometricalPredicates

# Arguments
# fd - funkcja odległości

struct Scaler
    scale::Any
    transx::Any
    transy::Any

    function Scaler(bbox::Array{Int,2})
        _scale = (abs(1 / (bbox[1] - bbox[2])))
        new(
            _scale,
            (1 + (0 - bbox[1] * _scale)),
            (1 + (0 - bbox[1, 2] * _scale)),
        )
    end
end

function scaled_point(unscaled_point::Point2D, scaler::Scaler)
    return Point(
        getx(unscaled_point) * scaler.scale + scaler.transx,
        gety(unscaled_point) * scaler.scale + scaler.transy,
    )
end

function unscaled_point(scaled_point::Point2D, scaler::Scaler)
    return Point(
        getx(scaled_point) / scaler.scale - scaler.transx,
        gety(scaled_point) / scaler.scale - scaler.transy,
    )
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
    h1 = calculateh1(h0)
    v1, v2 = extractbox(bbox, h0, h1)
    x, y = meshgrid(v1, v2)
    shiftevenrows!(x, h0)
    pb = [Point(row[1], row[2]) for row in eachrow([x[:] y[:]]) if fd(row) < geps]
    r0 = [1 ./ fh(getx(point), gety(point)) .^ 2 for point in pb]
    r0_max = maximum(r0 ./ maximum(r0))
    pold = Inf
    scaler = Scaler(bbox)
    scale = abs(1 / (bbox[1] - bbox[2]))
    transx = 1 + (0 - bbox[1] * scale)
    transy = 1 + (0 - bbox[1, 2] * scale)
    pb = [
        #Point(getx(point) * scale + transx, gety(point) * scale + transy)
        scaled_point(point, scaler)
        for point in pb if (rand(Float64, size(pb)))[1] < r0_max
    ]
    tesselation = DelaunayTessellation()
    push!(tesselation, pb)
    centers = Array{Point2D,1}()
    for triangle in tesselation
        center = centroid(Primitive(geta(triangle), getb(triangle), getc(triangle)))
        unscaled_p = unscaled_point(Point(getx(center), gety(center)), scaler)
        d = fd([getx(unscaled_p) gety(unscaled_p)])
        if d > -geps
            push!(centers, center)
        end
    end
    centers
    hbars = Array{Float64,1}()
    barvec = Array{Array{Float64,1},1}()
    for edge in delaunayedges(tesselation)
        b = unscaled_point(getb(edge), scaler)
        a = unscaled_point(geta(edge), scaler)
        # Obliczamy jako wartość fh() w połowie kazdego bara
        # push!(hbars, sqrt(length2(Line(a, b))))
        push!(
            barvec,
            [
              getx(b) - getx(a)
              gety(b) - gety(a)
            ],
        )
    end
    barvec
    L = [sqrt(sum(v_sum.^2)) for v_sum in barvec]
end

function dc(p)
    return sqrt(sum(p .^ 2)) - 1
end
