# compute momentum part

"""
X represents either below or above, i.e, < or > in the paper
esX, the array to energy point in region X, assuming equal probability
αX, the energy where nkX=0.5
βX, controls the fluctuation, βX is bigger when the local interaction is stronger
weightX, the weight of region X, in, X as <, weightX is nασ, X as >, weightX is 1-nασ
"""
function cal_nX(esX,αX,βX,weightX)
    nkX=[ cal_nk(αX,βX,esX_) for esX_ in esX]
    nX=mean(nkX)*weightX
end

"""
find the αX
# n_mean should be within 0.0 and 1.0
# or nX∈[0,weightX]
# nX=constraint_nX(nX,weightX)
αX_guess is computed by assuming esX has only one effective data point
"""
function solve_αX_from_nX(esX,nX,βX,weightX)
    n_mean=nX/weightX
    e_mean=mean(esX)
    αX_guess=cal_alpha(e_mean,βX,n_mean)
    αX=find_zero(αX->cal_nX(esX,αX,βX,weightX)-nX,αX_guess,Order0())
end

"""
this compute all necessary quantities for a given region
nX,AX,KX,nkX
nX : density (has been weighted)
AX : charge fluctuation (has been weighted)
KX : kinetic energy (has been weighted)
nkX: momentum density distribution
"""
function cal_nX_AX_KX_nkX(esX,αX,βX,weightX)
    nkX=[ cal_nk(αX,βX,esX_) for esX_ in esX]
    nX=mean(nkX)*weightX
    AX=mean(sqrt.(nkX.*(1.0 .- nkX)))*weightX
    KX=mean(esX.*nkX)*weightX
    nX,AX,KX,nkX
end


"""
compute momentum information
Abelow,Aabove,Kbelow,Kabove,[αbelow,αabove],[nkbelow,nkabove]

to break Aασ to Aασ_below and Aασ_above, so we can use the same form as the differential case. See ...v3.jl for details
nασ_,Δασ_,βασ_,eασ_=nασ[i],Δασ[i],βασ[i],eασ[i]
"""
function cal_Abelow_Aabove_Kbelow_Kabove_αασ_nk_(nασ_,Δασ_,βασ_,eασ_)
    βbelow,βabove=βασ_
    ebelow,eabove=eασ_
    nbelow=nασ_-Δασ_
    nabove=Δασ_
    αbelow=solve_αX_from_nX(ebelow,nbelow,βbelow,nασ_)
    αabove=solve_αX_from_nX(eabove,nabove,βabove,1-nασ_)
    nbelowcheck,Abelow,Kbelow,nkbelow=cal_nX_AX_KX_nkX(ebelow,αbelow,βbelow,nασ_)
    nabovecheck,Aabove,Kabove,nkabove=cal_nX_AX_KX_nkX(eabove,αabove,βabove,1-nασ_)
    Abelow,Aabove,Kbelow,Kabove,[αbelow,αabove],[nkbelow,nkabove]
end

# there are some function to compute derivatives, check them later. see scheme1xx

