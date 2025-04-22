@testitem "Aqua quality assurance tests" begin
    using Aqua: Aqua
    Aqua.test_all(
        PostNewtonian;
        ambiguities=false,
        unbound_args=(broken=true,),
        #persistent_tasks=(broken=true,),
    )
end
@testitem "ExplicitImports tests" begin
    using ExplicitImports: ExplicitImports
    @test ExplicitImports.check_no_implicit_imports(PostNewtonian) === nothing
    @test ExplicitImports.check_all_explicit_imports_via_owners(PostNewtonian) === nothing
    @test_broken ExplicitImports.check_all_explicit_imports_are_public(PostNewtonian) ===
        nothing
    @test ExplicitImports.check_no_stale_explicit_imports(
        PostNewtonian;
        ignore=(
            # The following generally need to be ignored because they are imported for use
            # in code that is created by a macro; as such, these are not explicitly used in
            # the code until after the macro runs (which comes after ExplicitImports looks
            # at the code).
            :M,
            :R,
            :X₁,
            :X₂,
            :q,
            :δ,
            :μ,
            :ν,
            :χ⃗₁,
            :χ⃗₂,
            :ℳ,
            :ln,
            :ln2,
            :ln3,
            :ln5,
            :order_index,
            :γₑ,
            :ζ3,
            :PNExpansionParameter,
            :type_converter,
            :MVector,
            :SVector,
            :Node,
        ),
    ) === nothing
    @test ExplicitImports.check_all_qualified_accesses_via_owners(PostNewtonian) === nothing
    @test_broken ExplicitImports.check_all_qualified_accesses_are_public(PostNewtonian) ===
        nothing
    @test ExplicitImports.check_no_self_qualified_accesses(
        PostNewtonian; ignore=(:Ω, :v)
    ) === nothing
end
