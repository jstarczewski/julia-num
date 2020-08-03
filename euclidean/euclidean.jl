function gcd(a::Int64, b::Int64)
    while b > 0
        c = a % b
        a = b
        b = c
    end
    return a
end
