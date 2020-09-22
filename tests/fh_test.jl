include("../mesh/mesh_clear.jl")

@testset "Simple values always same result" begin
    expected_result = 0.2
    @test fh(2, 3) == expected_result
    @test fh(0.2, 3) == expected_result
    @test fh(0.3, -0.1) == expected_result
end

@testset "Exception on wrong type of argument" begin
    @test_throws MethodError fh(1 + 1im, 3)
    @test_throws MethodError fh(1 + 2im, 3 + 3im)
    @test_throws MethodError fh(3 + 3im, 1 + 1im)
end

