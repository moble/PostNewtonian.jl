# Using this package from Python

While not quite as natural as calling Python code from Julia, it is very easy
to call Julia code from Python.  The process is essentially the same as using
any other Python package, other than installing Julia itself and any
dependencies *within* Julia that you may need (both of which are much easier
than similar tasks in Python).

!!! warning
    This package uses Unicode *internally* and in many of the examples found
    in this documentation.  However, every effort is made to ensure that
    *users* of this package can always use plain ASCII.  For example, if
    keyword arguments are accepted as Unicode, ASCII equivalents are usually
    also accepted.  Some function *names* are similarly in Unicode, but have
    ASCII equivalents.  See the various functions' documentation for
    acceptable replacements.

    It can be dangerous to use Unicode in Python in particular.  
    Python only accepts a small subset of Unicode — so that `M₁` for example
    is not valid input.  And what it does accept is automatically ["NFKC
    normalized"](https://en.wikipedia.org/wiki/Unicode_equivalence#Normal_forms).
    For example, variable names `Mₐ`, `Ma`, and `Mᵃ` are all treated
    identically by Python.  To illustrate this, consider the following
    python code:
    ```python
    >>> Mₐ = 1
    >>> Mᵃ = 2
    >>> Ma = 3
    >>> (Mₐ, Mᵃ, Ma)
    (3, 3, 3)
    ```
    We might have expected three different values `(1, 2, 3)` in the output,
    but Python never even sees the variable names as different strings; it
    interprets these expressions as setting, resetting, and resetting again
    the value of `Ma`.
    
    If you find an example where ASCII substitutions are not possible,
    please file a [bug
    report](https://github.com/moble/PostNewtonian.jl/issues/new).

Note that the Julia packages are installed uniquely to your python environment —
preferably a conda env.  For example, if you use two different conda envs to
call into Julia, you'll need to install the Julia packages for each env.  This
has the great advantage of allowing you to use different packages or versions in
each of the different environments.

1. Install Julia

   There are several options here:

   * If you are already use [conda](https://conda.io/) to manage your python
     installation (which is generally the recommended approach in python), and
     you are running one of the [platforms supported by
     conda-forge](https://anaconda.org/conda-forge/julia/files) (as of this
     writing, only `linux-64` or `osx-64` — not the newer ARM-based Macs, or any
     Windows machines), you can just add `julia` to the list of packages to
     install in step 2.
   * The [`juliaup`](https://github.com/JuliaLang/juliaup) installer is the most
     widely used option, and provides a simple method to update.
   * Otherwise, the official method is to just download a binary from the
     [official download page](https://julialang.org/downloads/).

   (It is also usually *very* easy to build Julia from source, but this should
   almost never be necessary.)

   Whichever method you choose, make sure that the `julia` executable is on your
   `PATH`.
   
2. Optionally, create a conda[^1] env just for this task
   ```bash
   conda create -c conda-forge -n julia_pn python numpy matplotlib
   conda activate julia_pn
   ```
   Add whatever other packages you use to that first line.

3. Install `juliacall` and `PostNewtonian`
   ```bash
   python -m pip install juliacall
   python -c 'from juliacall import Main as jl; jl.seval("""using Pkg; Pkg.add("PostNewtonian")""")'
   ```
   (Yes, you should use `pip` from *inside* a conda env.)  This will take a few
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

   # Declare some parameters
   delta = 0.0
   chi1 = np.array([0.1, 0.2, 0.3])
   chi2 = np.array([-0.2, 0.1, -0.3])
   Omega_i = 0.01

   # Call into Julia to run some function
   w = jl.GWFrames.PNWaveform("TaylorT1", delta, chi1, chi2, Omega_i)

   # Plot the magnitudes of all the modes as functions of time
   plt.semilogy(w.t, np.abs(w.data))
   ```
   The second-to-last line above uses the [`GWFrames.PNWaveform`](@ref) function
   from this package, which is meant to emulate the original syntax from the
   `GWFrames` package.  The resulting `w` will have various fields, like `t`,
   `data`, and `frame`, similar to those attached to `WaveformModes` objects in
   the `scri` and `sxs` packages.

In general, you can now call any Julia function by prepending `jl.` to the call
you would make in Julia.  As a fallback, you can evaluate actual Julia code in
the Julia session using `jl.seval("<Julia code goes here>")`.  This returns
whatever the Julia code would return.  A simple example is `x =
jl.seval("1+2")`.  See the [documentation for `juliacall`
here](https://github.com/cjdoris/PythonCall.jl#readme) for more details.

Typically, the main stumbling block is converting Python lists to Julia
`Vector`s when calling Julia functions.  Frequently, Julia code will have
difficulty if you try to pass a Python `list`, because `list`s do not have any
specific type that Julia can understand.  Instead, you should always convert a
list to a numpy array with `np.asarray`.  It is still possible that numpy will
not understand the type of the list, and you'll still get an error from Julia;
in this case you need to figure out [the `dtype` to tell numpy
about.](https://numpy.org/doc/stable/reference/generated/numpy.asarray.html).

Fortunately, conversion back to Python objects is more automatic.  In
particular, if Julia returns a `Vector`, `Matrix`, or more generally shaped
`Array`, you can usually just use that quantity in calls to Python functions.
If you really need a numpy array, the returned object will have a `.to_numpy()`
method.

Of course, it is *much* simpler to call Python code from Julia, so if you find
yourself using a lot of Julia code, you may want to consider flipping your
approach.


[^1]: As general advice, you should run `conda install -y mamba -n base -c
      conda-forge`, and then just use the command `mamba` wherever you would
      have used `conda`; `mamba` is a complete drop-in replacement, but is much
      faster because it's written in C instead of python.  For example, `mamba
      create -n julia_pn python numpy matplotlib` will typically run faster than
      the command given here.  This becomes a huge advantage when the env has
      lots of dependencies.
