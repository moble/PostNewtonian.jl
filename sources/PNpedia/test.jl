"""

We need to parse Mathematica expressions in Julia, using SymPy's parsing capabilities.
Unfortunately, at the moment, there is a bug in SymPy

    https://github.com/sympy/sympy/issues/27868

that prevents it from parsing, e.g., `Sqrt[2]*σ`, which does occur in these expressions.
There is a PR at

    https://github.com/sympy/sympy/pull/27876

that claims to fix this, but it has not been merged yet.

Also, David Trestini has evidently used `FullForm` to output the Mathematica to TXT files,
which translates the Greek letters to their full names — e.g., `ν` becomes `\[Nu]` — which
sympy cannot parse.  Specifically, `FullForm` encodes its output as "PrintableASCII".

There are some packages that may be useful for translating quite generally (in descending
order of how up-to-date I would guess they are):

    * `pygments-mathematica`
    * `mathics-scanner`
    * `FoxySheep2`

I'm wondering if it's worthwhile to to write a translator by getting all possible identifier
letters in Python, using

    import sys
    [c for c in map(chr, range(sys.maxunicode+1)) if ("_"+c).isidentifier()]

Then, we can take that into Mathematica, run each through `FullForm`, and create a JSON file
or something that we can read in with Python to replace the identifiers.  There are 139,463
such identifiers — some of which may not output at all with `FullForm`, and some of which
(not many) will be output exactly as input, so don't need to be included.

Note that Python's spec for identifiers lives at

    https://docs.python.org/3/reference/lexical_analysis.html#identifiers

A few useful pieces of information I glean from that:

> Within the ASCII range (U+0001..U+007F), the valid characters for identifiers include the
> uppercase and lowercase letters A through Z, the underscore _ and, except for the first
> character, the digits 0 through 9.

> All identifiers are converted into the normal form NFKC while parsing; comparison of
> identifiers is based on NFKC.

I suppose it's possible that some of these things that map to the same character in NFKC
could come up, and be encoded differently in Mathematica.  In that case, the output of the
"translate-to-python" function would be able to output different characters that map to the
same character in NFKC, but that's just a problem with how python works, so I don't know
that we can do anything about it.  Maybe warn about it?


One problem is unicode's combining characters, some of which

"""




using Pkg
Pkg.activate(@__DIR__)

import SymPyPythonCall
import PythonCall

greek_replacements = (
    raw"\[Nu]" => "ν",
    raw"\[Delta]" => "δ",
    raw"\[Tau]" => "τ",
    raw"\[Chi]" => "χ",
    raw"\[Sigma]" => "σ",
    raw"\[Kappa]" => "κ",
    raw"\[Lambda]" => "λ",
)



const parse_mathematica = PythonCall.pyimport("sympy.parsing.mathematica" => "parse_mathematica")

s = replace(
    raw"((1594323*I)/4480)*Sqrt[Pi/4199]*x^(9/2)*\[Nu]*(\[Delta] - 6*\[Delta]*\[Nu] + 10*\[Delta]*\[Nu]^2 - 4*\[Delta]*\[Nu]^3)",
    greek_replacements...
)

ex = parse_mathematica(s)#, Dict("EllipticE[x]"=>"elliptic_e(x)"))



const elliptic_e = PythonCall.pyimport("sympy.functions.special.elliptic_integrals" => "elliptic_e")

s = raw"(Sqrt[Pi])/(2*EllipticE[m])"
s = raw"((1594323*I)/4480)*Sqrt[Pi/4199]*x^(9/2)*\[Nu]*(\[Delta] - 6*\[Delta]*\[Nu] + 10*\[Delta]*\[Nu]^2 - 4*\[Delta]*\[Nu]^3)"
ex = parse_mathematica(s)#, Dict("EllipticE[x]"=>"elliptic_e(x)"))

ex.args[2].args[0].args
