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


## Installation
 

### 0. Optionally pre-install Julia

If you'll want to use Julia from outside of Python, install Julia
first with [the `juliaup`
installer](https://github.com/JuliaLang/juliaup?tab=readme-ov-file#juliaup---julia-version-manager).

### 1. Install the `sxs` python package

While not *strictly necessary*, the `sxs` python package is very
useful for providing a general interface to waveforms, and will
automatically install the PostNewtonian package and its dependencies
(and even Julia, if necessary).

If you use `conda`, just run something like[^1]
```bash
conda install -c conda-forge sxs numba::numba numba::llvmlite
```
If you prefer `pip`, use
```bash
python -m pip install sxs
```

!!! danger "Avoid segfaults on macOS"

    If using conda/mamba on macOS, you *must* install the correct
    versions of `numba` and `llvmlite`, or [you will get
    segfaults.](https://github.com/numba/numba/issues/7857#issuecomment-1082246028)
    Specifically, you *must* include the `numba::` prefix to
    select the correct anaconda channel, as shown in the
    command given above.  A more permanent solution is to run the
    following commands:
    ```
    conda config --set channel_priority strict
    conda config --add channels conda-forge --prepend channels numba
    ```
    This will automatically choose the correct versions of numba and
    llvmlite for your system, without the need to explicitly use the
    `numba::` prefix ever again.

### 2. Install Julia and `PostNewtonian.jl`

The following line will take a few minutes to install Julia (if
necessary), then download and compile all the necessary Julia
packages.

```bash
python -c 'from sxs import julia'
```


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

# Call `sxs` wrapper for Julia functions
w = PNWaveform(M1, M2, chi1, chi2, Omega_i)

# Plot the magnitudes of all the modes as functions of time
plt.semilogy(w.t, np.abs(w.data))
```
The [`sxs.julia.PNWaveform`
function](https://github.com/sxs-collaboration/sxs/blob/e6aa63695fdb1a2f97cfb54e04dbbd5453142cd3/sxs/julia/__init__.py#L17-L86)
is a simple wrapper that calls [`orbital_evolution`](@ref) and
[`inertial_waveform`](@ref), then wraps the result in an
[`sxs.WaveformModes`
object](https://sxs.readthedocs.io/en/stable/api/waveforms/#waveformmodes-class),
with all the usual methods and properties, as well as several
additional fields relevant to post-Newtonian waveforms:

- `M1` (array): The primary mass as a function of time.
- `M2` (array): The secondary mass as a function of time.
- `chi1` (array): The primary spin as a function of time.
- `chi2` (array): The secondary spin as a function of time.
- `frame` (array): The quaternionic frame as a function of time.
- `v` (array): The orbital velocity as a function of time.
- `orbital_phase` (array): The orbital phase as a function of
    time.

The input arguments for
[`sxs.julia.PNWaveform`](https://github.com/sxs-collaboration/sxs/blob/e6aa63695fdb1a2f97cfb54e04dbbd5453142cd3/sxs/julia/__init__.py#L17-L86)
are mostly the same as for [`orbital_evolution`](@ref), except that
the keyword argument `inertial` will switch whether the returned
waveform will be in the inertial frame (the default of `True`) or the
coorbital frame (`False`), and the keyword arguments `ell_min`,
`ell_max`, and `waveform_pn_order` will be intercepted and passed to
[`inertial_waveform`](@ref) or [`coorbital_waveform`](@ref)
correspondingly.


## Full Julia interface from Python

In general, you can now call any function from the Julia
`PostNewtonian` package by running
```python
from sxs.julia import PostNewtonian
```
and then calling the function just as you would if you had run `import
PostNewtonian` in Julia.  That is, any of the functions given in the
rest of this documentation should be available in Python, just as they
would be in Julia.  (Again, note that all function names and arguments
should have ASCII equivalents to use from Python.)

You can also evaluate *arbitrary Julia code* by prefixing
`PostNewtonian.` to the command, or as a fallback write Julia code as
a string and pass it to `PostNewtonian.seval("<Julia code goes
here>")`.  Simple examples include
```python
>>> PostNewtonian.println("Hello from Julia!")
Hello from Julia!
>>> x = PostNewtonian.seval("1+2")
>>> x
3
```
See the [documentation for `juliacall`
here](https://github.com/cjdoris/PythonCall.jl#readme) for more
details.

Whenever they correspond naturally to built-in Python types — like
`int`, `float`, `str`, etc. — the returned objects will be converted
to such Python objects.  Otherwise, the returned objects will be
wrappers around Julia `Vector`s, `Matrix`es, or more generally shaped
`Array`s.  Usually, you will be able to pass these objects directly to
Python functions.  If you really need a `numpy` array, you can use the
`to_numpy()` method that the wrappers will support.

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
