# Using this package from Python

While not quite as natural as calling Python code from Julia, it is
very easy to call Julia code from Python.  The process is essentially
the same as using any other Python package, other than installing
Julia itself and any dependencies *within* Julia that you may need
(both of which are much easier than similar tasks in Python).

!!! warning "Be careful of using unicode in Python!"

    Every effort is made to ensure that *users* of this package can always use
    plain ASCII — even though this package uses Unicode *internally* and in many
    of the examples found in this documentation.  For example, if keyword
    arguments are accepted as Unicode, ASCII equivalents are usually also
    accepted.  Some function *names* are similarly in Unicode, but have ASCII
    equivalents.  See the various functions' documentation for acceptable
    replacements.

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
    >>> print((Mₐ, Mᵃ, Ma))
    (3, 3, 3)
    ```
    We might have expected three different values `(1, 2, 3)` in the output,
    but Python never even sees the variable names as different strings; it
    interprets these expressions as setting, resetting, and resetting again
    the value of `Ma`.
    
    If you find an example where ASCII substitutions are not possible,
    please file a [bug
    report](https://github.com/moble/PostNewtonian.jl/issues/new).


## Getting started

Installation is pretty simple.  There are two optional steps first,
but generally just one required step.  

!!! tip "TL;DR"
    Basically, you just need to install the `sxs` package, and then
    run the python command `from sxs import julia`.


### Optional step 1: Install Julia

If a Julia installation is not already available, the made below will
automatically install it for you *in the python environment you are
using*.  If you have multiple python environments, this means that
there will be multiple copies of Julia installed, one in each
environment.  Installing Julia beforehand will avoid this duplication.

The officially recommended method is to use [the `juliaup`
installer](https://github.com/JuliaLang/juliaup?tab=readme-ov-file#juliaup---julia-version-manager),
which also provides an easy way to update Julia itself, manage
multiple versions of Julia, and update it.

Note that the `julia` executable must be on your `PATH` for the next
steps to work properly.  You may need to restart your shell or
terminal after installing Julia to make sure.


### Optional step 2: Create a conda env

Create a conda[^1] env just for this task
```bash
conda create -c conda-forge -n julia_pn python numpy matplotlib
conda activate julia_pn
```
Add whatever other packages you use to that first line.


### Required step: Install the `sxs` package and use it to install Julia

While the `sxs` package is not *strictly* required, it is very useful
for providing a general interface to waveforms, and will automatically
install the `PostNewtonian` package and its dependencies.  Simply run

```bash
python -m pip install sxs
python -c 'from sxs import julia'
```
The first line installs `sxs` itself, but not Julia; the second line
will take a few minutes to install Julia then download and compile all
the necessary packages.


## Testing the installation
   
Start up a python session and run something like this:
```python
# Any python imports you need go here
import numpy as np
import matplotlib.pyplot as plt

# Start the Julia session
from sxs.julia import PNWaveform

# Declare the essential parameters
M1 = 0.5
M2 = 0.5
chi1 = [0.1, 0.2, 0.3]
chi2 = [-0.2, 0.1, -0.3]
Omega_i = 0.01

# Call into Julia to run some function
w = PNWaveform(M1, M2, chi1, chi2, Omega_i)

# Plot the magnitudes of all the modes as functions of time
plt.semilogy(w.t, np.abs(w.data))
```
The `w` object returned here will be an `sxs.WaveformModes` object,
with all the usual methods and properties.

---

In general, you can now call any function from the Julia
`PostNewtonian` package by running
```python
from sxs.julia import PostNewtonian
```
and then calling the function just as you would if you had run `import
PostNewtonian` in Julia.  As a fallback, you can evaluate actual Julia
code in the Julia session using `PostNewtonian.seval("<Julia code goes
here>")`.  This returns whatever the Julia code would return.  A
simple example is `x = PostNewtonian.seval("1+2")`.  See the
[documentation for `juliacall`
here](https://github.com/cjdoris/PythonCall.jl#readme) for more
details.

Sometimes, the returned objects will be wrappers around Julia
`Vector`s, `Matrix`es, or more generally shaped `Array`s.  Usually,
you will be able to pass these objects directly to Python functions.
If you really need a `numpy` array, you can use the `to_numpy()`
method that the wrappers will support.

Of course, it is *much* simpler to call Python code from Julia, so if
you find yourself using a lot of Julia code, you may want to consider
flipping your approach.


[^1]: As general advice, you should run `conda install -y mamba -n
      base -c conda-forge`, and then just use the command `mamba`
      wherever you would have used `conda`; `mamba` is a complete
      drop-in replacement, but is much faster because it's written in
      C instead of python.  For example, `mamba create -n julia_pn
      python numpy matplotlib` will typically run faster than the
      command given here.  This becomes a huge advantage when the env
      has lots of dependencies.
