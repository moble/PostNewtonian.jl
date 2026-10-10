@testitem "differentiability" begin
    using ForwardDiff
    using PostNewtonian
    using PostNewtonian: symbol_index, symbols
    using Random
    Random.seed!(1234)

    # First, test that the gradient and Hessian of ℰ′ have zeros whenever a spin, Rotor, or
    # Φ component is involved and the PN order is ≤1

    pn = rand(PostNewtonian.BBH, PNOrder=2//2)
    function c1(u)
        pnsystem = PostNewtonian.BBH(u; PNOrder=PostNewtonian.pn_order(pn))
        return PostNewtonian.𝓔′(pnsystem)
    end
    ∇c = ForwardDiff.gradient(c1, pn.state)
    Hc = ForwardDiff.hessian(c1, pn.state)

    pn_type = typeof(pn)
    nonzero_indices = [symbol_index(pn_type, s) for s ∈ [:M₁, :M₂, :v]]
    zero_indices = [
        i for i ∈ eachindex(symbols(pn_type)) if i ∉ nonzero_indices
    ]

    @test c1(pn.state) == PostNewtonian.𝓔′(pn)
    @test all(y->y≠0, ∇c[nonzero_indices])
    @test all(y->y==0, ∇c[zero_indices])
    @test all(ij->Hc[ij...]≠0, Iterators.product(nonzero_indices, nonzero_indices))
    @test all(ij->Hc[ij...]==0, Iterators.product(zero_indices, zero_indices))

    # Next, test that the Hessian of ℰ′ has zeros whenever a spin or Φ component is involved
    # and the PN order is 3//2.  This is the order at which linear-in-spin terms appear, but
    # no quadratic-in-spin terms.
    pn = rand(PostNewtonian.BBH, PNOrder=3//2)
    function c2(u)
        pnsystem = PostNewtonian.BBH(u; PNOrder=PostNewtonian.pn_order(pn))
        return PostNewtonian.𝓔′(pnsystem)
    end
    ∇c = ForwardDiff.gradient(c2, pn.state)
    Hc = ForwardDiff.hessian(c2, pn.state)

    pn_type = typeof(pn)
    zero_indices = [
        symbol_index(pn_type, :χ⃗₁ˣ):symbol_index(pn_type, :χ⃗₁ᶻ);
        symbol_index(pn_type, :χ⃗₂ˣ):symbol_index(pn_type, :χ⃗₂ᶻ);
        symbol_index(pn_type, :Φ)
    ]
    nonzero_indices = [
        i for i ∈ eachindex(symbols(pn_type)) if i ∉ zero_indices
    ]

    @test c2(pn.state) == PostNewtonian.𝓔′(pn)
    @test all(y->y≠0, ∇c[1:(end - 1)])
    @test all(ij->Hc[ij...]≠0, Iterators.product(nonzero_indices, nonzero_indices))
    @test all(ij->Hc[ij...]==0, Iterators.product(zero_indices, zero_indices))

    # Now test that the gradient and Hessian are entirely nonzero when the PN order is >1,
    # except for the very last component, which is the phase, which shouldn't affect
    # anything
    pn = rand(PostNewtonian.BBH, PNOrder=4//2)
    function c3(u)
        pnsystem = PostNewtonian.BBH(u; PNOrder=PostNewtonian.pn_order(pn))
        return PostNewtonian.𝓔′(pnsystem)
    end
    ∇c = ForwardDiff.gradient(c3, pn.state)
    Hc = ForwardDiff.hessian(c3, pn.state)
    @test c3(pn.state) == PostNewtonian.𝓔′(pn)
    @test all(y->y≠0, ∇c[1:(end - 1)])
    @test ∇c[end] == 0
    @test all(y->y≠0, Hc[1:(end - 1), 1:(end - 1)])
    @test all(y->y==0, Hc[1:end, end])
    @test all(y->y==0, Hc[end, 1:end])

    # Finally, check that the full-order energy and waveforms are also differentiable
    pn = rand(PostNewtonian.BBH; v=0.5)
    function c4(u)
        pnsystem = PostNewtonian.BBH(u; PNOrder=PostNewtonian.pn_order(pn))
        return PostNewtonian.𝓔′(pnsystem)
    end
    ∇c = ForwardDiff.gradient(c4, pn.state)
    Hc = ForwardDiff.hessian(c4, pn.state)
    @test c4(pn.state) == PostNewtonian.𝓔′(pn)
    @test all(y->y≠0, ∇c[1:(end - 1)])
    @test ∇c[end] == 0
    @test all(y->y≠0, Hc[1:(end - 1), 1:(end - 1)])
    @test all(y->y==0, Hc[1:end, end])
    @test all(y->y==0, Hc[end, 1:end])
    function c5(u)
        pnsystem = PostNewtonian.BBH(u; PNOrder=PostNewtonian.pn_order(pn))
        inspiral = PostNewtonian.orbital_evolution(pnsystem; vₑ=0.6)
        h = PostNewtonian.inertial_waveform(inspiral)
        return sum(abs2, h)
    end
    ∇c = ForwardDiff.gradient(c5, pn.state)
    Hc = ForwardDiff.hessian(c5, pn.state)
    @test c5(pn.state) > 0
    @test all(y->y≠0, ∇c[1:(end - 1)])
    @test ∇c[end] == 0
    @test all(y->y≠0, Hc[1:(end - 1), 1:(end - 1)])
    @test all(y->y==0, Hc[1:end, end])
    @test all(y->y==0, Hc[end, 1:end])
end
