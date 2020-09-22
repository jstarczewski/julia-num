using Plots
x = 1:10; y = rand(10); # These are the plotting data
plot(x, y);

println("HELLO")


function f(z)
    a = sum(z.^2,dims=2)
    return a
end

println(f([1 2; 3 4]))

l = [1 2; 3 4]
println(l[1,1])
println(l[2,1])
println(l[1,2])

for a in 1:10
end

bbox = [0 0; 2 2]
h0 = 1
x = bbox[1,1]:h0:bbox[2,1]
y = bbox[1,2]:h0:bbox[2,2]
println(x)
println(y)

# Funkcja meshgrid znana z Matlaba jest niedostępna w Julii
# Tworzoymy jej podstawową implementację
# Obliczamy długość podanych w parametrach funkcji wektorów
# Reprezentują one wymiarny naszej siatki
# Nastepnie funkcja repeat tworzy wektor skladajacy sie z
# Elemenntów pierwotnego wektora, ale kazdy element zostal powtórzon????????????y m ilość razy
# Co finalnie daje nam długi wektor wspolrzednych x punktów naszej siatki
# Funkcja reshape łamie wektor do macierzy o wymiarach m na na
# Ta sama prodecura powtórzona jest dla punktów y tylko w tym przypadku
# Funkcja repeat powtarza cały wektor y n razy, na koniec łamie go do macierzy
# Zwracana jest para macierzy gx, gy


function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where T
    m, n = length(vy), length(vx)

    gx = reshape(repeat(vx, inner = m, outer = 1), m, n)
    gy = reshape(repeat(vy, inner = 1, outer = n), m, n)

    return gx, gy
end

println(repeat([1 2; 3 4], inner=(2,1), outer=(1, 3)))
println(repeat([1 2; 3 4], inner=(2,3)))
#println(meshgrid(x))

p = 