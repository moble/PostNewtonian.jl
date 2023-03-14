# Automatic differentiation

While symbolic differentiation can be useful in many important scenarios, it is
not very helpful for providing derivatives of general evolved waveforms, because
waveforms require integration of ODEs.  However, Julia is quite capable of
automatically differentiating even through an ODE integration — both gradients
and Hessians — with respect to any or all of the initial conditions.

There are many automatic differentiation packages in Julia.  Zygote will most
likely be the best in the not-too-distant future, but at the time of writing,
the simplest approach is to use
[`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl).  A simple
wrapper function that takes only the arguments to be differentiated may be
needed.  The returned quantity will be a vector of waveforms corresponding to
the derivative in the waveform at each instant with respect to the desired
parameters.

More broadly, it can also be helpful to differentiate a *function* of a
waveform.  For example, if we are trying to minimize the difference between a
waveform and a PN waveform, we may have a cost function that takes (some or all
of) the PN parameters integrates the ``L^2`` norm of the difference between
them.  This cost function should be easily differentiable as well.
