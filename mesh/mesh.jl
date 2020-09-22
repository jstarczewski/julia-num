using Gadfly
using LinearAlgebra
using VoronoiDelaunay
using GeometricalPredicates

function distmesh(fd, fh, h0, bbox, pfix) end

# Co jest ten broadcasting aka .+ w 12

function fh(x::T, y::T) where {T}
    return 0.2
end

function generate(fd, fh, bbox, h0, pfix)
    ttol = 0.1
    geps = 0.001 * h0
    h1 = h0 * (sqrt(3) / 2)
    v1 = bbox[1, 1]:h0:bbox[2, 1]
    v2 = bbox[1, 2]:h1:bbox[2, 2]
    x, y = meshgrid(v1, v2)
    x[2:2:end, :] = x[2:2:end, :] .+ h0 / 2
    p = [x[:] y[:]]
    z = [vcat(a) for a in eachrow(p) if fd(a) < geps]
    f = transpose(reshape(vcat(z...), 2, length(z)))
    r0 = 1 ./ fh(f[:, 1], f[:, 2]) .^ 2
    r0_max = r0 ./ max(r0)
    l = [vcat(b) for b in eachrow(f) if (rand(Float64, size(f)))[1] < r0_max]
    g = transpose(reshape(vcat(z...), 2, length(z)))
    g
    N = size(g, 1)
    N
    println(N)

    pold = Inf
    scale = abs(1 / (bbox[1] - bbox[2]))
    transx = 1 + (0 - bbox[1] * scale)
    transy = 1 + (0 - bbox[1, 2] * scale)

    maximum(sqrt.(sum((g .- pold) .^ 2, dims = 2) / h0)) > ttol
    tess = DelaunayTessellation()
    aa = map(el -> Point(el[1] * scale + transx, el[2] * scale + transy), eachrow(g))
    aa
    println(aa)

    push!(tess, aa)
    zda = Array{Point2D,1}()
    for tri in tess
        push!(zda, centroid(Primitive(geta(tri), getb(tri), getc(tri))))
    end
    zza = filter(cent -> dc([getx(cent) gety(cent)]) > -geps, zda)
    zdda = Array{Float64,1}()
    barvec = Array{Array{Float64,1},1}()
    for edge in delaunayedges(tess)
        push!(zdda, sqrt(length2(Line(geta(edge), getb(edge)))) / scale)
        push!(
            barvec,
            [
                (getx(getb(edge)) - getx(geta(edge))) / scale,
                ((gety(getb(edge)) - gety(geta(edge)))) / scale,
            ],
        )
    end
    hbars = 0.2
    ss = 0.2 * 1.2 * sqrt(sum(zdda .^ 2) / (0.2^2))
    ttt = Array{Float64,1}()
    for edge in delaunayedges(tess)
        push!(ttt, ss - sqrt(length2(Line(geta(edge), getb(edge)))) / scale)
    end
    rttt = reshape(ttt, :, 1)
    mzdda = reshape(zdda, :, 1)
    init_Fvec = rttt ./ mzdda
    #    Fvec = Fvec.*barvec
    l
    #    println(sort(with_f[:]))

    #mzzda
    #    width = max_coord - min_coord
    #    a = Point2D[Point(min_coord + rand() * width, min_coord + rand() * width) for i in 1:100]
    #    push!(tess, a)
    #    x, y = getplotxy(delaunayedges(tess))
    #    Gadfly.plot(x=x, y=y, Geom.path)
end



function evaluate(fd, p, geps)
    for x in eachrow(p)
        if (fd(x) < geps)

        end
    end
end

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where {T}
    m, n = length(vy), length(vx)

    gx = reshape(repeat(vx, inner = m, outer = 1), m, n)
    gy = reshape(repeat(vy, inner = 1, outer = n), m, n)

    return gx, gy
end

function feval(f, x...) end

function generateNodes()
    x = [1, 2, 3, 4]
    y = [1, 1, 2, 3]
    plot(x, y)
end

function dc(p)
    return sqrt(sum(p .^ 2)) - 1
end

function dcircle(p, xc, yc, r)
    return sqrt.((p[:, 1] .- xc) .^ 2 .+ (p[:, 2] .- yc) .^ 2) .- r
end
