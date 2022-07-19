# GWFrames compatibility layer

[`GWFrames`](https://github.com/moble/GWFrames) is an old python package that
was too hard to upgrade to python 3, and therefore became too hard to maintain.
Many of its capabilities have been superseded by the
[`scri`](https://github.com/moble/scri) and
[`sxs`](https://github.com/sxs-collaboration/sxs/) packages, but not the
post-Newtonian capabilities.  This package now supersedes the latter.  For
convenience, the `PostNewtonian.GWFrames` submodule provides a simple function
to mimic the arguments used with the original `GWFrames` package to obtain a PN
waveform (with a couple extra arguments).  We won't bother to provide a return
type that can fully mimic the object returned by the original `GWFrames`
package, though it should be similar to the ones used by `scri` and `sxs`.

```@autodocs
Modules = [PostNewtonian]
Pages   = ["gwframes.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```
