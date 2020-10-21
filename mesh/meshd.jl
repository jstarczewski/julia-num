using Deldir
using Gadfly

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


function dc(p)
    #return -0.3 + abs(0.7 - sqrt(sum(p .^ 2)))
    return sqrt(sum(p .^ 2)) - 1
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
    p = [x[:] y[:]]
    p = [vcat(row) for row in eachrow(p) if fd(row) < -geps]
    r0 = [1 ./ fh(row...) .^ 2 for row in p]
    r0_max = maximum(r0 ./ maximum(r0))
    pold = Inf
    scaler = Scaler(bbox)
    p = [vcat(scale_p(row, scaler)) for row in p if (rand(Float64, size(p)))[1] < r0_max]
    p = transpose(reshape(vcat(p...), 2, length(p)))
    del, vor, summ = deldir(p[:, 1], p[:, 2])
    Dx, Dy = edges(del)
    interior_points = Array{Array{Float64}, 1}()
    p = [Dx[:], Dy[:]]
    generators = Dict{Int64, Array{Int,1}}()
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
           i = findall(in(generators[e]),generators[k])
           common = generators[k][i]
           for c in common
              tri = sort([c, k, e])
              if !(tri in triangles)
                   push!(triangles, tri)
               end
           end
       end
    end
    println(triangles)
    Gadfly.plot(x=Dx, y=Dy, Geom.path, Coord.cartesian(fixed=true))
end
