module PointScaler

using GeometricalPredicates

export Scaler, scale_p, unscale_p, scaled_point, unscaled_point

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
        ((getx(unscaled_point) * scaler.scale) + scaler.transx),
        ((gety(unscaled_point) * scaler.scale) + scaler.transy),
    )
end

function unscaled_point(scaled_point::Point2D, scaler::Scaler)
    return Point(
        (getx(scaled_point) - scaler.transx) / scaler.scale,
        (gety(scaled_point) - scaler.transy) / scaler.scale,
    )
end

end