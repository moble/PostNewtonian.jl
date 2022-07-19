# Using this package from python


1. Install Julia
   
   The [`juliaup`](https://github.com/JuliaLang/juliaup) installer has pretty
   broad support from some of the main Julia developers, and seems likely to be
   the main installation method in the not-too-distant future.
   
   Otherwise, the official method is to just download a binary from the
   [official download page](https://julialang.org/downloads/).
   
   Make sure that the `julia` binary is on your `PATH`.
   
2. Create a conda[^1] env just for this task
   ```bash
   conda create -n julia_pn python numpy matplotlib
   conda activate julia_pn
   ```
   Add whatever other packages you use to that first line.

3. Install `juliacall` and build `PostNewtonian`
   ```bash
   python -m pip install juliacall
   python -c 'from juliacall import Main as jl; jl.seval("""using Pkg; Pkg.add("PostNewtonian")""")'
   ```
   This will take a few minutes to compile all the necessary packages.

4. Test the installation.  Start up a python session and run something like
   this:
   ```python
   import numpy as np
   import matplotlib.pyplot as plt
   from juliacall import Main as jl, convert as jlconvert
   
   VectorF64 = jl.seval("Vector{Float64}")
   jl.seval("using PostNewtonian")
   delta = 0.0
   chi1 = jlconvert(VectorF64, [0.1, 0.2, 0.3])
   chi2 = jlconvert(VectorF64, [-0.2, 0.1, -0.3])
   Omega_i = 0.01
   w = jl.GWFrames.PNWaveform("TaylorT1", delta, chi1, chi2, Omega_i)
   ```
   The resulting `w` will have various fields, like `t`, `data`, and `frame`,
   similar to those attached to `WaveformModes` objects in the `scri` and `sxs`
   packages.

The last line above uses the [`PNWaveform`](@ref) function from this package,
which is meant to emulate the original syntax from the `GWFrames` package.


[^1]: In general, it's better to `conda install -y mamba` in your "base" conda
      env, and then just use `mamba` wherever you would have used `conda`;
      `mamba` is a complete drop-in replacement, but is much faster because
      it's written in C instead of python.
