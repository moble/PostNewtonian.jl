@testset verbose=true "Aqua quality assurance tests" begin
    Aqua.test_all(PostNewtonian; ambiguities=false)
end
@testset verbose=true "ExplicitImports tests" begin
    @test ExplicitImports.check_no_implicit_imports(PostNewtonian) === nothing
    @test ExplicitImports.check_all_explicit_imports_via_owners(PostNewtonian) === nothing
    #@test ExplicitImports.check_all_explicit_imports_are_public(PostNewtonian) === nothing
    @test ExplicitImports.check_no_stale_explicit_imports(
        PostNewtonian;
        ignore=(
            :M, :R, :X₁, :X₂, :q, :δ, :μ, :ν, :χ⃗₁, :χ⃗₂, :ℳ,
            :ln, :ln2, :ln3, :ln5, :order_index, :γₑ, :ζ3
        ),
    ) === nothing
    @test ExplicitImports.check_all_qualified_accesses_via_owners(PostNewtonian) === nothing
    #@test ExplicitImports.check_all_qualified_accesses_are_public(PostNewtonian) === nothing
    @test ExplicitImports.check_no_self_qualified_accesses(
        PostNewtonian;
        ignore=(:Ω, :v,)
    ) === nothing
end
