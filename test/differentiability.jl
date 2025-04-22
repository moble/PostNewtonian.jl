@testitem "differentiability" begin
    using ForwardDiff
    using PostNewtonian

    using Random
    Random.seed!(1234)

    # First, test that the gradient and Hessian of ℰ′ have zeros whenever a spin, Rotor, or
    # Φ component is involved and the PN order is ≤1
    nonzero_indices = [PostNewtonian.M₁index, PostNewtonian.M₂index, PostNewtonian.vindex]
    zero_indices = [
        i for i ∈ eachindex(PostNewtonian.pnsystem_symbols) if i ∉ nonzero_indices
    ]
    pn = rand(PostNewtonian.BBH, PNOrder=2//2)
    function c(u)
        pnsystem = PostNewtonian.BBH(u; PNOrder=PostNewtonian.pn_order(pn))
        return PostNewtonian.𝓔′(pnsystem)
    end
    ∇c = ForwardDiff.gradient(c, pn.state)
    Hc = ForwardDiff.hessian(c, pn.state)
    @test c(pn.state) == PostNewtonian.𝓔′(pn)
    @test all(y->y≠0, ∇c[nonzero_indices])
    @test all(y->y==0, ∇c[zero_indices])
    @test all(ij->Hc[ij...]≠0, Iterators.product(nonzero_indices, nonzero_indices))
    @test all(ij->Hc[ij...]==0, Iterators.product(zero_indices, zero_indices))

    # Next, test that the Hessian of ℰ′ has zeros whenever a spin or Φ component is involved
    # and the PN order is 3//2.  This is the order at which linear-in-spin terms appear, but
    # no quadratic-in-spin terms.
    zero_indices = [
        PostNewtonian.χ⃗₁indices;
        PostNewtonian.χ⃗₂indices;
        PostNewtonian.Φindex
    ]
    nonzero_indices = [
        i for i ∈ eachindex(PostNewtonian.pnsystem_symbols) if i ∉ zero_indices
    ]
    pn = rand(PostNewtonian.BBH, PNOrder=3//2)
    function c(u)
        pnsystem = PostNewtonian.BBH(u; PNOrder=PostNewtonian.pn_order(pn))
        return PostNewtonian.𝓔′(pnsystem)
    end
    ∇c = ForwardDiff.gradient(c, pn.state)
    Hc = ForwardDiff.hessian(c, pn.state)
    @test c(pn.state) == PostNewtonian.𝓔′(pn)
    @test all(y->y≠0, ∇c[1:(end - 1)])
    @test all(ij->Hc[ij...]≠0, Iterators.product(nonzero_indices, nonzero_indices))
    @test all(ij->Hc[ij...]==0, Iterators.product(zero_indices, zero_indices))

    # Now test that the gradient and Hessian are entirely nonzero when the PN order is >1,
    # except for the very last component, which is the phase, which shouldn't affect
    # anything
    pn = rand(PostNewtonian.BBH, PNOrder=4//2)
    function c(u)
        pnsystem = PostNewtonian.BBH(u; PNOrder=PostNewtonian.pn_order(pn))
        return PostNewtonian.𝓔′(pnsystem)
    end
    ∇c = ForwardDiff.gradient(c, pn.state)
    Hc = ForwardDiff.hessian(c, pn.state)
    @test c(pn.state) == PostNewtonian.𝓔′(pn)
    @test all(y->y≠0, ∇c[1:(end - 1)])
    @test ∇c[end] == 0
    @test all(y->y≠0, Hc[1:(end - 1), 1:(end - 1)])
    @test all(y->y==0, Hc[1:end, end])
    @test all(y->y==0, Hc[end, 1:end])

    # Finally, check that the full-order energy and waveforms are also differentiable
    pn = rand(PostNewtonian.BBH; v=0.5)
    function c(u)
        pnsystem = PostNewtonian.BBH(u; PNOrder=PostNewtonian.pn_order(pn))
        return PostNewtonian.𝓔′(pnsystem)
    end
    ∇c = ForwardDiff.gradient(c, pn.state)
    Hc = ForwardDiff.hessian(c, pn.state)
    @test c(pn.state) == PostNewtonian.𝓔′(pn)
    @test all(y->y≠0, ∇c[1:(end - 1)])
    @test ∇c[end] == 0
    @test all(y->y≠0, Hc[1:(end - 1), 1:(end - 1)])
    @test all(y->y==0, Hc[1:end, end])
    @test all(y->y==0, Hc[end, 1:end])
    function c(u)
        pnsystem = PostNewtonian.BBH(u; PNOrder=PostNewtonian.pn_order(pn))
        inspiral = PostNewtonian.orbital_evolution(pnsystem; vₑ=0.6)
        h = PostNewtonian.inertial_waveform(inspiral)
        return sum(abs2, h)
    end
    ∇c = ForwardDiff.gradient(c, pn.state)
    Hc = ForwardDiff.hessian(c, pn.state)
    @test c(pn.state) > 0
    @test all(y->y≠0, ∇c[1:(end - 1)])
    @test ∇c[end] == 0
    @test all(y->y≠0, Hc[1:(end - 1), 1:(end - 1)])
    @test all(y->y==0, Hc[1:end, end])
    @test all(y->y==0, Hc[end, 1:end])
end
