# This should be an internal class that just allows us to combine ODE solution objects
struct CombinedInterpolationData <: SciMLBase.AbstractDiffEqInterpolation
    sol₋
    sol₊
    tᵢ
end
function (interp::CombinedInterpolationData)(tvals::AbstractFloat, idxs, deriv, p, continuity::Symbol=:left)
    if tvals ≥ interp.tᵢ
        interp.sol₊.interp(tvals, idxs, deriv, p, continuity)
    else
        interp.sol₊.interp(tvals, idxs, deriv, p, continuity)
    end
end
function (interp::CombinedInterpolationData)(val, tvals::AbstractFloat, idxs, deriv, p, continuity::Symbol=:left)
    if tvals ≥ interp.tᵢ
        interp.sol₊.interp(val, tvals, idxs, deriv, p, continuity)
    else
        interp.sol₊.interp(val, tvals, idxs, deriv, p, continuity)
    end
end
function (interp::CombinedInterpolationData)(tvals, idxs, deriv, p, continuity::Symbol=:left)
    val = if idxs isa Integer
        Vector{eltype(interp.sol₊)}(undef, length(tvals))
    elseif isnothing(idxs)
        [Vector{eltype(interp.sol₊)}(undef, size(interp.sol₊, 1)) for _ in eachindex(tvals)]
    else
        [Vector{eltype(interp.sol₊)}(undef, length(idxs)) for _ in eachindex(tvals)]
    end
    interp(val, tvals, idxs, deriv, p, continuity)
end
function (interp::CombinedInterpolationData)(val, tvals, idxs, deriv, p, continuity::Symbol=:left)
    i₊ = tvals .≥ interp.tᵢ
    i₋ = tvals .< interp.tᵢ
    if idxs isa Integer
        for i in eachindex(tvals)
            if tvals[i] ≥ interp.tᵢ
                val[i] = interp.sol₊.interp(tvals[i], idxs, deriv, p, continuity)
            else
                val[i] = interp.sol₋.interp(tvals[i], idxs, deriv, p, continuity)
            end
        end
    else
        for i in eachindex(tvals)
            if tvals[i] ≥ interp.tᵢ
                val[i] .= interp.sol₊.interp(tvals[i], idxs, deriv, p, continuity)
            else
                val[i] .= interp.sol₋.interp(tvals[i], idxs, deriv, p, continuity)
            end
        end
    end
    syms = interp.sol₊.prob.f.syms
    indepsym = interp.sol₊.prob.f.indepsym
    observed = interp.sol₊.prob.f.observed
    p = interp.sol₊.prob.p
    RecursiveArrayTools.DiffEqArray(val, tvals, syms, indepsym, observed, p)
end


"""
    combine_solutions(sol₋, sol₊)

Combine ODESolutions

This function is internal to this package.  It is not entirely general, but allows us to
combine the backwards- and forwards-in-time solutions of the PN orbital-evolution ODE
equations into a single `ODESolution` object that should behave just as if it were the
result of `solve`.  In particular, indexing, interpolation, and iterations should behave
exactly as [described in the `DifferentialEquations`
docs](https://diffeq.sciml.ai/stable/basics/solution/).

"""
function combine_solutions(sol₋, sol₊)
    alg = sol₊.alg
    t = [reverse(sol₋.t[2:end]); sol₊.t]
    u = [reverse(sol₋.u[2:end]); sol₊.u]
    retcode = sol₊.retcode  # Could be something more clever; maybe the worse retcode?
    problem = ODEProblem(sol₊.prob.f, u[1], (t[1], t[end]), sol₊.prob.p)
    if sol₊.dense
        interp = CombinedInterpolationData(sol₋, sol₊, sol₊.t[1])
        DiffEqBase.build_solution(problem, alg, t, u, dense=true, retcode=retcode, interp=interp)
    else
        DiffEqBase.build_solution(problem, alg, t, u, dense=false, retcode=retcode)
    end
end
