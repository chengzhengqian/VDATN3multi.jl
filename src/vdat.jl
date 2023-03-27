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

# compute obserbles and some useful utility function
include("./utils.jl")

# the code is generated from mathematica
include("./include_gene_code.jl")

# utilies for discretizating band structure
include("./load_band.jl")

# compute local matrix, (we have move the functions to construct w using fluctuation into a different file)
include("./local.jl")

# using the fluctuation approach to construct w
# the problem is that there is a boundary,
# try the exponetial form, also, add symmetry
include("./w_fluc.jl")
# there is another approach, use exponential form
# previously, we only use the uniform transform to adjust the total energy, this approach could also be used with exponeital form, thus provide a reliable way to constraint the density for the multiband case
# !! we remove the uniform transform part, this is mainly for the free density case.  
include("./w_free.jl")

# momentum part; we should remove the derivatives part, check the original version for how to compute derivatives
include("./momentum.jl")


# generate two particle interaction
include("./interation.jl")

# data structure for model
include("./model.jl")

include("./compute.jl")

