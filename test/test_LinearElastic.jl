
@testset "test_LinearElastic" begin
    
    E = 210e3
    ν = 0.3

    # Test for 2D elastic material
    mat_2D = LinearElastic{2, Float64}(E, ν)
    ε_2D = SymTensorValue{2, Float64}(1e-3, 0, 0)
    σ_2D = mat_2D.C_elas ⊙ ε_2D
    @test σ_2D isa SymTensorValue{2, Float64}

    state0 = initial_material_state(mat_2D)
    σ_new, C_elas, state_new = material_response(mat_2D, ε_2D, state0)
    @test σ_new isa SymTensorValue{2, Float64}
    @test C_elas isa SymFourthOrderTensorValue{2, Float64}
    @test state_new isa LinearElasticState{2, Float64}

    # Test for 3D elastic material
    mat_3D = LinearElastic{3, Float64}(E, ν)
    ε_3D = SymTensorValue{3, Float64}(1e-3, 0, 0, 0, 0, 0)
    σ_3D = mat_3D.C_elas ⊙ ε_3D
    @test σ_3D isa SymTensorValue{3, Float64}

    state0 = initial_material_state(mat_3D)
    σ_new, C_elas, state_new = material_response(mat_3D, ε_3D, state0)
    @test σ_new isa SymTensorValue{3, Float64}
    @test C_elas isa SymFourthOrderTensorValue{3, Float64}
    @test state_new isa LinearElasticState{3, Float64}

end

