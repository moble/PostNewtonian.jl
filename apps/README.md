From thee directory containing this file, run the following commands
to create a simple Julia application, compile it into a standalone
executable using PackageCompiler.jl, and run it:

```bash
> julia -q --project=.

julia> using PackageCompiler

julia> create_app("ZeroEccParamsFromPN", "ZeroEccParamsFromPNApp")
PackageCompiler: bundled libraries:
...
✔ [00m:35s] PackageCompiler: compiling nonincremental system image

julia> exit()

> ZeroEccParamsFromPNApp/bin/ZeroEccParamsFromPN
```
