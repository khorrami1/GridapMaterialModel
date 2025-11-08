
@testset "test_wrappers" begin
    
    E = 210e3
    ν = 0.3
    
    C_elas_2D = elastic_tangent(2, E, ν)
    @test C_elas_2D isa SymFourthOrderTensorValue{2, Float64}
    
    C_elas_3D = elastic_tangent(3, E, ν)
    @test C_elas_3D isa SymFourthOrderTensorValue{3, Float64}
    
    ε_2D = SymTensorValue{2, Float64}(1e-3, 0, 0)
    @test C_elas_2D ⊙ ε_2D isa SymTensorValue{2, Float64}

    ε_3D = SymTensorValue{3, Float64}(1e-3, 0, 0, 0, 0, 0)
    @test C_elas_3D ⊙ ε_3D isa SymTensorValue{3, Float64}

end