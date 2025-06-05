To compile this app, you need Julia installed.  Then, from the
directory containing this file, run the following commands to create a
simple Julia application, compile it into a standalone executable
using PackageCompiler.jl, and run it:

```bash
> julia --threads=auto make.jl

[lots of output...]
[might take ~five minutes...]

> ZeroEccParamsFromPNApp/bin/ZeroEccParamsFromPN --q=4.3 --chiA=0.1,0.2,0.3 --chiB=0.3,0.2,0.1 --Omega0=0.015
```

In particular, be sure to run `ZeroEccParamsFromPN` at least once,
because it may have to install some pre-built libraries the first time
you run it.  At this point, you can move the `ZeroEccParamsFromPNApp`
directory to a different location, and/or make it read-only, so that
other users can run it without having to install Julia or any
dependencies.

Note that the app will be quite large (~750 MB) when created with
Julia <1.12, but this is expected to shrink significantly in future
versions of Julia.
