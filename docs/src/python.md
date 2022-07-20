# Using this package from Python

While not quite as natural as calling Python code from Julia, it is very easy
to call Julia code from Python.  The process is essentially the same as using
any other Python package, other than installing Julia itself and any
dependencies *within* Julia that you may need (both of which are much easier
than similar tasks in Python).

Note that the Julia packages are installed uniquely to your python environment
— preferably a conda env.  For example, if you use two different conda envs to
call into Julia, you'll need to install the Julia packages for each env.  This
has the great advantage of allowing you to use different packages or versions
in each of the different environments.

1. Install Julia
   
   The [`juliaup`](https://github.com/JuliaLang/juliaup) installer has pretty
   broad support from some of the main Julia developers, and seems likely to be
   the main installation method in the not-too-distant future.
   
   Otherwise, the official method is to just download a binary from the
   [official download page](https://julialang.org/downloads/).
   
   Make sure that the `julia` executable is on your `PATH`.
   
2. Optionally, create a conda[^1] env just for this task
   ```bash
   conda create -n julia_pn python numpy matplotlib
   conda activate julia_pn
   ```
   Add whatever other packages you use to that first line.

3. Install `juliacall` and install `PostNewtonian`
   ```bash
   python -m pip install juliacall
   python -c 'from juliacall import Main as jl; jl.seval("""using Pkg; Pkg.add("PostNewtonian")""")'
   ```
   (Yes, you should use `pip` from inside a conda env.)  This will take a few
   minutes to compile all the necessary packages in Julia.

4. Test the installation
   
   Start up a python session and run something like this:
   ```python
   # Any python imports you need go here
   import numpy as np
   import matplotlib.pyplot as plt
   
   # Start the Julia session
   from juliacall import Main as jl
   
   # Import `PostNewtonian` in the Julia session
   jl.seval("using PostNewtonian")

   # Convenience function for converting to a Julia `Vector`
   Vector = lambda a: jl.convert(jl.Vector, a)

   # Declare some parameters
   delta = 0.0
   chi1 = Vector([0.1, 0.2, 0.3])
   chi2 = Vector([-0.2, 0.1, -0.3])
   Omega_i = 0.01
   
   # Call into Julia to run some function
   w = jl.GWFrames.PNWaveform("TaylorT1", delta, chi1, chi2, Omega_i)
   ```
   The last line above uses the [`GWFrames.PNWaveform`](@ref) function from
   this package, which is meant to emulate the original syntax from the
   `GWFrames` package.  The resulting `w` will have various fields, like `t`,
   `data`, and `frame`, similar to those attached to `WaveformModes` objects in
   the `scri` and `sxs` packages.

In general, you can now call any Julia function by prepending `jl.` to the call
you would make in Julia.  As a fallback, you can evaluate actual Julia code in
the Julia session using `jl.seval("<Julia code goes here>")`.  This returns
whatever the Julia code would return.  A simple example is `x =
jl.seval("1+2")`.  See the [documentation for `juliacall`
here](https://github.com/cjdoris/PythonCall.jl#readme) for more details.

The main stumbling block is converting Python lists and numpy arrays to Julia
arrays when calling Julia functions.  In the code above, we see an example of
how to do this conveniently: by defining a Python function `Vector` that
converts Python lists (or 1-D numpy arrays) to Julia `Vector`s.  Note that this
will try to automatically infer the type of the list — `float` in the case
shown above.  If you want to be more specific, you can use Julia's standard
parametric types, replacing braces `{}` with brackets `[]`.  For example, to
specify that the conversion should be to the Julia type `Vector{Float64}`, we
could define
```python
VectorF64 = lambda a: jl.convert(jl.Vector[jl.Float64], a)
```
Note that `Vector([0, 0, 0])` will result in a `Vector{Int64}`, whereas
`VectorF64([0, 0, 0])` will result in a `Vector{Float64}`.

Fortunately, conversion back to Python objects is more automatic.  In
particular, if Julia returns a `Vector`, you can usually just use that quantity
in calls to Python functions.  If you really need a numpy array, the returned
object will have a `.to_numpy()` method.

Of course, it is *much* simpler to call Python code from Julia, so if you find
yourself using a lot of Julia code, you may want to consider flipping your
approach.


[^1]: As general advice, you should `conda install -y mamba` in your "base"
      conda env, and then just use the command `mamba` wherever you would have
      used `conda`; `mamba` is a complete drop-in replacement, but is much
      faster because it's written in C instead of python.  For example, `mamba
      create -n julia_pn python numpy matplotlib` will typically run faster
      than the command given here.  This becomes a huge advantage when the env
      has lots of dependencies.
