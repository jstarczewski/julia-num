function triangles(del::DataFrame, summ::DataFrame)
    generators = emptygenerators(summ)
    generatorstogenerated!(generators, del)
    triangles = buildtriangles(generators)
    return map(
        t -> [
            Point(summ[t[1], 1], summ[t[1], 2]),
            Point(summ[t[2], 1], summ[t[2], 2]),
            Point(summ[t[3], 1], summ[t[3], 2]),
        ],
        triangles,
    )
end

function emptygenerators(summ::DataFrame)::Dict{Int64,Array{Int,1}}
    generators = Dict{Int64,Array{Int,1}}()
    for (index, row) in enumerate(eachrow(summ))
        generators[index] = []
    end
    return generators
end

function generatorstogenerated!(generators::Dict{Int64,Array{Int,1}}, del::DataFrame)
    for row in eachrow(del)
        generators[row[5]] = append!(generators[row[5]], row[6])
        generators[row[6]] = append!(generators[row[6]], row[5])
    end
end

function buildtriangles(generators::Dict{Int64,Array{Int,1}})::Array{Array{Int,1},1}
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
    return triangles
end
