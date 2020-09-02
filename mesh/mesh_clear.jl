using Test

# Arguments
# fd - funkcja odległości

function fh(x::T, y::T)::Real where {T<:Real}
    @info "Value is" s = 0.2
    return 0.2
end

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where {T<:Real}
    m, n = length(vy), length(vx)
    @info "Length m = " m
    @info "Length n = " n
    gx = reshape(repeat(vx, inner = m, outer = 1), m, n)
    gy = reshape(repeat(vy, inner = 1, outer = n), m, n)
    @info "Grid x" gx
    @info "Grid y" gy
    return gx, gy
end

#bbox = [0 0; 1 1]

function extractvectorsfromboundry(bbox::Array{Real,2}, h0::Real, h1::Real)
    v1 = bbox[1, 1]:h0:bbox[2, 1]
    v2 = bbox[1, 2]:h1:bbox[2, 2]
    return v1, v2
end

function calculateh1(h0::Real)
    return h0 * (sqrt(3) / 2)
end

function generate(fd, fh, bbox, h0, pfix)

    h1 = calculateh1(h0)
    v1, v2 = extractvectorsfromboundry(bbox, h0, h1)
    x, y = meshgrid(v1, v2)
end
