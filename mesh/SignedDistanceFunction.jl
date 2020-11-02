
function dc(p)
    return -0.3 + abs(0.7 - sqrt(sum(p .^ 2)))
end

function dca(p)
    return sqrt(sum(p .^ 2)) - 1
end
