
@testset "test_J2Plasticity" begin
    
    E = 210e3
    ν = 0.3
    σ0 = 210.0
    H = 10e-1

    # Test for 2D materials
    material_2D = J2Plasticity{2, Float64}(E, ν, σ0, H)
    ε = SymTensorValue{2, Float64}(5e-3, 0., 0.)
    state0 = initial_material_state(material_2D)
    σ_new, C_ep, state_new = material_response(material_2D, ε, state0, nothing, nothing, nothing)
    @test σ_new isa SymTensorValue{2, Float64}
    @test C_ep isa SymFourthOrderTensorValue{2, Float64}
    @test state_new isa J2PlasticityState{2, Float64}

    # Test for 3D material
    material_3D = J2Plasticity{3, Float64}(E, ν, σ0, H)
    @test material_3D.C_elas isa SymFourthOrderTensorValue{3, Float64}
    ε = SymTensorValue{3, Float64}(5e-3, 0., 0., 0., 0., 0.)
    state0 = initial_material_state(material_3D)
    σ_new, C_ep, state_new = material_response(material_3D, ε, state0, nothing, nothing, nothing)
    @test σ_new isa SymTensorValue{3, Float64}
    @test C_ep isa SymFourthOrderTensorValue{3, Float64}
    @test state_new isa J2PlasticityState{3, Float64}

end

##############################################################
# Test one element 

# using Plots

# function uniaxialTest(loadingRange, Δε)
#     m = J2Plasticity{3, Float64}(;E=210e3, ν=0.3, σ0=200.0, H=10.0e3)
#     #cache = get_cache(m)
#     state = initial_material_state(m)
#     σ_equivalent_all = Float64[]
#     σ_all = SymTensorValue{3, Float64}[]
#     state_all = J2PlasticityState[]
#     ε_all = SymTensorValue{3, Float64}[]
#     push!(σ_equivalent_all, 0.0)
#     push!(σ_all, zero(SymTensorValue{3, Float64}))
#     push!(ε_all, zero(SymTensorValue{3, Float64}))
#     push!(state_all, state)
#     ε = zero(SymTensorValue{3,Float64})
#     for e11 in loadingRange
#         ε += Δε
#         σ, ∂σ∂ε, state = material_response(m, ε, state, nothing, nothing, nothing)
#         push!(σ_equivalent_all, sqrt(1.5*dev(σ)⊙dev(σ)))
#         push!(σ_all, σ)
#         push!(state_all, state)
#         push!(ε_all, ε)
#     end
#     return ε_all, σ_all, σ_equivalent_all, state_all
# end

# loadingRange = range(0.0, 0.02, 201)
# Δε = SymTensorValue{3, Float64}( 0.0, loadingRange.step.hi, 0.0, 0.0, 0.0, 0.0)
# ε_all, σ_all, σ_equivalent_all, state_all = uniaxialTest(loadingRange, Δε)
# # p = plot([e[1,2] for e in ε_all], [s[1,2] for s in σ_all])
# # p = plot([e[1,2] for e in ε_all], [s for s in σ_equivalent_all])
# loadingRange = range(0.0, 0.02, 201)
# Δε = SymTensorValue{3,Float64}(loadingRange.step.hi, 0.0, 0.0, -0.5*loadingRange.step.hi, 0.0, -0.5*loadingRange.step.hi)
# ε_all, σ_all, σ_equivalent_all, state_all = uniaxialTest(loadingRange, Δε)
# # p = plot([e[1,1] for e in ε_all], [s[1,1] for s in σ_all])
# # p = plot([e[1,1] for e in ε_all], [s for s in σ_equivalent_all])

# # End of one element test
##############################################################


