# Apps and scripts

Apps and scripts are meant to be used from the command line, rather
than a running Julia session.  The main difference is that an app is a
standalone executable, while a script requires a working installation
of Julia and a project with `PostNewtonian` installed.  Apps can also
be moved around on a system and used by other users without requiring
them to install Julia or any dependencies, and will be slightly
(usually a few seconds) faster to run because everything will have
been compiled ahead of time.  The scripts are geared more towards
developers or users of the `PostNewtonian` package itself.

## ZeroEccParamsFromPN

As of this writing, the only app supplied by this package is
`ZeroEccParamsFromPN`, and is also available as a script.  The reason
you would use it as a script is if you want to use it with the current
version of `PostNewtonian` — if you are developing the package or
something like that.

To compile the app, follow the directions in
`PostNewtonian/apps/README.md`.  Run either the app or the script with
the `--help` option to see the available options.

## Other scripts

Primarily as conveniences for developers, two other scripts are
included: `docs.jl` and `test.jl`.

The `docs.jl` script is used to generate and serve the documentation
for the package.  This can be helpful because it will automatically
rebuild the documentation if you make changes to the documentation
markdown pages — although you need to restart the script to see
changes made to docstrings.

The `test.jl` script is used to run the tests for the package, and
processes coverage reports.
