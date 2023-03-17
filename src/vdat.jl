# using the G12, u, β scheme for N=3. Allow general density-density interaction
# see VDAT.jl, and the code on ginsburg for more details

using LinearAlgebra
using Statistics
using Combinatorics
using Roots
using NLsolve
using Optim
using Dierckx
using DelimitedFiles

include("./utils.jl")

# the code is generated from mathematica
include("./include_gene_code.jl")

# utilies for discretizating band structure
include("./load_band.jl")

include("./local.jl")
