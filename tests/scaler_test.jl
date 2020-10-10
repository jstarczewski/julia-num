# include("../mesh/mesh_clear.jl")

struct Scaler
    scale::Any
    transx::Any
    transy::Any

    function Scaler(bbox::Array{Int,2})
        _scale = (abs(1 / (bbox[2] - bbox[1 ])))
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
        (getx(scaled_point) - scaler.transx) / scaler.scale,
        (gety(scaled_point) - scaler.transy) / scaler.scale,
    )
end

bbox = [-1 -1; 1 1]
@info "Boundy box = " bbox
@info "x1 = " bbox[1]
@info "y1 = " bbox[2]

_scale = (abs(1 / (bbox[1] - bbox[2])))

@info "Scale = " _scale