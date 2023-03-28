```@meta
CurrentModule = PostNewtonian
```

# PostNewtonian

This package computes orbital dynamics of and waveforms from binary black-hole
systems, in the [post-Newtonian
approximation](https://en.wikipedia.org/wiki/Post-Newtonian_expansion).
Currently, there is no support for eccentric systems, but general precessing
quasicircular systems are supported.

```@example
using PlotlyJS
p = PlotlyJS.Plot(sin.(1:0.1:10))
open("sin.html","w") do f
    PlotlyJS.PlotlyBase.to_html(f, p; include_plotlyjs="cdn", full_html=true)
end
```

```@raw html
<iframe src="sin.html" style="height:500px;width:100%;"></iframe>
```

```julia
# @example
using Plots
plotlyjs()  # Requires adding `PlotlyJS` to your project

savefig("sol.html")  # hide
```
```@raw html
<iframe src="sol.html" style="height:500px;width:100%;"></iframe>
```
