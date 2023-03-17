include("./gene_code/c_func.jl")
include("./gene_code/pmatw.jl")
include("./gene_code/g11matw.jl")
include("./gene_code/g11matwS.jl")
include("./gene_code/g12matw.jl")
include("./gene_code/g12matwS.jl")
include("./gene_code/g33matw.jl")
include("./gene_code/g33matwS.jl")
# just check
include("./gene_code/g13matw.jl")
include("./gene_code/g23matw.jl")
include("./gene_code/g31matw.jl")
include("./gene_code/g32matw.jl")
include("./gene_code/neffInn.jl")
include("./gene_code/gloc12NonInteractingInnloc.jl")
include("./gene_code/sloc11sloc12.jl")
# we also add the function for fixed point
include("./gene_code/rIns11s12.jl")
include("./gene_code/chalf.jl")
include("./gene_code/g012Incr.jl")
include("./gene_code/g011Incr.jl")
include("./gene_code/nk.jl")
# we add derivatives now
include("./gene_code/dnAdab.jl")
include("./gene_code/alpha.jl")
include("./gene_code/glocInnA.jl")
# we may need a function that only compute the remaining part
# i.e g12,g13,g23,g31,g32             
include("./gene_code/glocReducedInA.jl")
include("./gene_code/delta.jl")
include("./gene_code/g0IngG12.jl")
include("./gene_code/g0fullIngG12.jl")

# check for fixed point (examples of compute g012 from fixed Δ and S)
# nloc=0.4
# g012=0.3
# g12=0.2
# sloc1,sloc2=cal_sloc11sloc12(nloc,g012,g12)
# r=cal_rIns11s12(sloc1,sloc2)
# Δ=cal_delta(g12,sloc1,sloc2)
# c=cal_chalf(nloc,Δ,sloc1,sloc2) # this is wrong, use r are guments
# c=cal_chalf(nloc,Δ,r) # this is the correct version
# g012_check=cal_g012Incr(nloc,Δ,r,c)
# g011_check=cal_g011Incr(nloc,Δ,r,c)
# just A block of g0, 
# g012=0.4                        

# cal_pmatw(g012)
# cal_g11matw(g012)
# cal_g11matwS(g012)
# cal_g12matw(g012)
# cal_g12matwS(g012)
# # we need full entries of g0
# g012,g013,g023,g031,g032,g033=rand(6)
# cal_g33matw(g012,g013,g023,g031,g032,g033)
# cal_g33matwS(g012,g013,g023,g031,g032,g033)
# #
# nloc=0.3
# cal_neffInn(nloc,g012)
# cal_gloc12NonInteractingInnloc(nloc,g012)
# nloc=0.5
# g12=0.1
# g012=0.5
# cal_sloc11sloc12(nloc,g012,g12)
# α=0.2
# β=0.3
# ϵ=0.4
# cal_nk(α,β,ϵ)
# the output order is {dnkda,dnkdb,dAkda,dAkdb}
# cal_dnAdab(α,β,ϵ)
# # args={epsilon,beta,nk};
# nk=0.999;
# nk=0.001;
# cal_alpha(ϵ,β,nk)
# nloc=0.3
# delta=0.1
# Abelow=Aabove=0.1
# sloc11=0.2
# sloc12=0.3
# # function cal_glocInnA(n,delta,Abl,Aab,sloc11,sloc12)
# # order of arguments
# # gloc=cal_glocInnA(nloc,delta,Abelow,Aabove,sloc11,sloc12)
# # glocReduced=cal_glocReducedInA(Abelow,Aabove,sloc11,sloc12)
# # [gloc[1,3],gloc[2,3],gloc[3,1],gloc[3,2]].-glocReduced
# # cal_delta(g12,sloc11,sloc12)
# g13=g23=0.2
# g31=g32=-0.2
# cal_g0IngG12(nloc,g012,g12,g13,g23,g31,g32)
