# GWFrames compatibility layer

[`GWFrames`](https://github.com/moble/GWFrames) is an old python
package that was too hard to upgrade to python 3, and therefore became
impossible to maintain.  Many of its capabilities have been superseded
by the [`scri`](https://github.com/moble/scri) and
[`sxs`](https://github.com/sxs-collaboration/sxs/) packages, but not
the post-Newtonian capabilities.  This package now supersedes the
latter.  For convenience, the `PostNewtonian.GWFrames` submodule
provides a simple function to mimic the arguments used with the
original `GWFrames` package to obtain a PN waveform (with a couple
extra arguments).  We won't bother to provide a return type that can
fully mimic the object returned by the original `GWFrames` package,
though it should be similar to the ones used by `scri` and `sxs`.

This interface can be used directly from Julia, but is more likely to
be wanted by python users.  To borrow from the example in [the
standard python interface documentation](@ref
Testing-the-installation):
```julia
# Any python imports you need go here
import numpy as np
import matplotlib.pyplot as plt

# Start the Julia session
from sxs.julia import GWFrames

# Declare the essential parameters
delta = 0.0
chi1 = [0.1, 0.2, 0.3]
chi2 = [-0.2, 0.1, -0.3]
Omega_i = 0.01

# Call into Julia to run some function
w = GWFrames.PNWaveform("TaylorT1", delta, chi1, chi2, Omega_i)

# Plot the magnitudes of all the modes as functions of time
plt.semilogy(w.t, np.abs(w.data))
```
This should produce precisely the same plot as the python example, but
using the `GWFrames` interface.  Note that the parameters given here
are required, but optional parameters are the same as in the original
`GWFrames` interface — which are listed in full in the docstring
below.

When called directly from Julia, the result will just be a
`NamedTuple` with fields described in the docstring below.  When
called from python as shown above, the resulting `w` will be an
[`sxs.WaveformModes`
object](https://sxs.readthedocs.io/en/stable/api/waveforms/#waveformmodes-class),
with all the usual methods and properties, as well as several fields
relevant to post-Newtonian waveforms:

- `M1` (array): The primary mass as a function of time.
- `M2` (array): The secondary mass as a function of time.
- `chi1` (array): The primary spin as a function of time.
- `chi2` (array): The secondary spin as a function of time.
- `frame` (array): The quaternionic frame as a function of time.
- `v` (array): The orbital velocity as a function of time.
- `orbital_phase` (array): The orbital phase as a function of
    time.


```@autodocs
Modules = [PostNewtonian.GWFrames]
Pages   = ["gwframes.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```
