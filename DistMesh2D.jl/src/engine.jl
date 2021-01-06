"""
    Funnkcja generujaca krawędzie jest punktem wspolnym obu silnikow.
    Generacja krawedzi wymaga generacji trojkatow, które w zaleznosci od
    użytego algorytmu są generowane w mniej lub bardziej złożone sposób
"""
function dd(
    p,
    iteration::Int64,
    fd::Function,
    scaler::Scaler,
    geps::Float64
)::Array{GeometricalPredicates.Line2D{GeometricalPredicates.Point2D},1}
    """
        Kod zewnętrznej biblioteki jest zamkniety tylko w tej funkcji
    """
    del, vor, summ = deldir(p[:, 1], p[:, 2])
    """
        W DelDir generacja trojkatow jest zaimplementowana przeze mnie na podstawie
        położenia "generatorów" (punktów generacji)
    """
    trigs = triangles(del, summ)
    edges =  [Line(Point(r[1], r[2]), Point(r[3], r[4])) for r in eachrow(del)]
    return validedges(trigs, edges, scaler, fd, geps)
end
#Obie funkcje zwracają taki sam typ danych, ale nie da sie go narzucić strukturze
function vd(
    p,
    iteration::Int64,
    fd::Function,
    scaler::Scaler,
    geps::Float64
)::Array{GeometricalPredicates.Line2D{GeometricalPredicates.Point2D},1}
    """
        Kod zewnętrznej biblioteki jest zamkniety tylko w tej funkcji
    """
    tess = DelaunayTessellation(convex = true)
    p = [Point2D(pp[1], pp[2]) for pp in eachrow(p)]
    push!(tess, p)
    """
        W VoronoiDelaunay.jl trójkaty sa od razu stworzone, wystarczy sie po nich
        przeiterowac i zmapowac do interesujacej nas struktury
    """
    trigs = triangles(tess)
    edges = Array{GeometricalPredicates.Line2D{GeometricalPredicates.Point2D},1}()
    # Nie wspiera iteratora
    for edge in delaunayedges(tess)
        push!(edges, Line(geta(edge), getb(edge)))
    end
    return validedges(trigs, edges, scaler, fd, geps)
end

# Abstrakcyjny wspolny typ, dla poprawnosci
abstract type Engine end

"""
    Silnik DelDir, zbudowany z odpowiedniego scalera, bo Biblioteka pracuje na zakresie
    bbox = [0 0; 1 1]
"""
struct DD <: Engine
    scaler::Scaler
    edges::Function

    function DD(bbox::Array{Float64, 2})
        _scaler = Scaler(bbox)
        # Przekazanie odpowiedniej funkcji generujacej
        new(Scaler(bbox), dd)
    end
end

#Silnik VoronoiDelauny, oparty o bbox = [1 1; 2 2]
struct VD <: Engine
    scaler::Scaler
    edges::Function

    function VD(bbox::Array{Float64, 2})
        # Przekazane przesuniecie do scalera
        new(Scaler(bbox, basetrans = 1.0), vd)
    end
end
