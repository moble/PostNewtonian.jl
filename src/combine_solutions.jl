# This should be an internal class that just allows us to combine ODE solution objects
struct CombinedInterpolationData <: SciMLBase.AbstractDiffEqInterpolation
    sol₋
    sol₊
    tᵢ
end
function (interp::CombinedInterpolationData)(tvals, idxs, deriv, p, continuity::Symbol=:left)
    val = similar(interp.sol₊.u, length(tvals))
    interp(val, tvals, idxs, deriv, p, continuity)
end
function (interp::CombinedInterpolationData)(val, tvals, idxs, deriv, p, continuity::Symbol=:left)
    i₊ = tvals .≥ interp.tᵢ
    i₋ = tvals .< interp.tᵢ
    if any(i₊)
        val₊ = interp.sol₊.interp(tvals[i₊], idxs, deriv, p, continuity)
        val[i₊] .= val₊.u
    end
    if any(i₋)
        val₋ = interp.sol₋.interp(tvals[i₋], idxs, deriv, p, continuity)
        val[i₋] .= val₋.u
    end
    RecursiveArrayTools.DiffEqArray(val, tvals, val₊.syms, val₊.indepsym, val₊.observed, val₊.p)
end


"""
    combine_solutions(sol₋, sol₊)

Combine ODESolutions

This function is internal to this package.  It is not entirely general, but
allows us to combine the backwards- and forwards-in-time solutions of the PN
inspiral ODE equations into a single `ODESolution` object that should behave
just as if it were the result of `solve`.  In particular, indexing,
interpolation, and iterations should behave exactly as [described in the
`DifferentialEquations` docs](https://diffeq.sciml.ai/stable/basics/solution/).

"""
function combine_solutions(sol₋, sol₊)
    alg = sol₊.alg
    t = [reverse(sol₋.t); sol₊.t]
    u = [reverse(sol₋.u); sol₊.u]
    retcode = sol₊.retcode  # Could be something more clever; maybe the worse retcode?
    problem = ODEProblem(sol₊.prob.f, u[1], (t[1], t[end]))
    if sol₊.dense
        interp = CombinedInterpolationData(sol₋, sol₊, sol₊.t[1])
        DiffEqBase.build_solution(problem, alg, t, u, dense=true, retcode=retcode, interp=interp)
    else
        DiffEqBase.build_solution(problem, alg, t, u, dense=false, retcode=retcode)
    end
end
