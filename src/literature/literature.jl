# These files define functions encapsulating common variables used throughout the
# literature.  Note that not every reference uses these conventions, which is why
# they must be explicitly imported to be used by `@pn_reference` modules.
include("common_variables/horizons.jl")
include("common_variables/mass_combinations.jl")
include("common_variables/orbital_elements.jl")
include("common_variables/spin_combinations.jl")
include("common_variables/tidal_coupling.jl")

# Now we include the `@pn_reference` modules themselves, named by the bibtex keys we use for
# them in ``../references.bib`.
#include("references/Einstein1918.jl")
