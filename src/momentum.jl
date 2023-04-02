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
# debug when nασ->1
1-nασ_-Δασ_
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


# in the previous way, we first  assume the eασX has only one point and estimate the αX, then solve αX. This works well for normal value, but for n->0 or n->1, this has singularities. So in these cases, we could simply assume that nk is flat (or batter approximation?)
"""
we merge solving αX with cal_nX_AX_KX_nkX
# new version
esX=[1,2,3]
nX=1e-8
βX=4.0
weightX=0.1+1e-8
cal_nX_AX_KX_nkX_with_nX(esX,nX,βX,weightX)
"""
function cal_nX_AX_KX_nkX_with_nX(esX,nX,βX,weightX)
    αX=solve_αX_from_nX(esX,nX,βX,weightX)
    cal_nX_AX_KX_nkX(esX,αX,βX,weightX)
end

"""
assumig nX_mean ->1, or 0, i.e, in this limit, we can safely assume the density distribution is flat and equals nX_mean
cal_nX_AX_KX_nkX_with_nX_mean(1.0,esX,weightX)
cal_nX_AX_KX_nkX_with_nX_mean(0.0,esX,weightX)
# also, this function
cal_nX_AX_KX_nkX_with_nX_mean(0.0,esX,0.0)
"""
function cal_nX_AX_KX_nkX_with_nX_mean(nX_mean,esX,weightX)
    nkX=[nX_mean for esX_ in esX]
    cal_nX_AX_KX_nkX_from_nkX(nkX,esX,weightX)
end

function cal_nX_AX_KX_nkX_from_nkX(nkX,esX,weightX)
    nX=mean(nkX)*weightX
    AX=mean(sqrt.(nkX.*(1.0 .- nkX)))*weightX
    KX=mean(esX.*nkX)*weightX
    nX,AX,KX,nkX
end

"""
this one consider some limiting behaverisou
cal_nX_AX_KX_nkX_safe(esX,0.0,1.0,0.5)
cal_nX_AX_KX_nkX_safe(esX,0.5,1.0,0.5)
cal_nX_AX_KX_nkX_safe(esX,0.499999,1.0,0.5)
cal_nX_AX_KX_nkX_safe(esX,0.0,1.0,0.5)
cal_nX_AX_KX_nkX_safe(esX,0.0001,1.0,0.2)
"""
function cal_nX_AX_KX_nkX_safe(esX,nX,βX,weightX)
    if(weightX==0.0)
        return cal_nX_AX_KX_nkX_with_nX_mean(0.5,esX,0.0)
    else
        nX_mean=nX/weightX
        if(nX_mean==0.0 || nX_mean==1.0)
            return cal_nX_AX_KX_nkX_with_nX_mean(nX_mean,esX,weightX)
        else
            return cal_nX_AX_KX_nkX_with_nX(esX,nX,βX,weightX)
        end        
    end    
end

"""
including the limiting case, we don't save αασ now, we can provide a function to compute it if the user wants.
nασ_,Δασ_,βασ_,eασ_=nασ[i],Δασ[i],βασ[i],eασ[i]
"""
function cal_Abelow_Aabove_Kbelow_Kabove_nk_safe(nασ_,Δασ_,βασ_,eασ_)
    βbelow,βabove=βασ_
    ebelow,eabove=eασ_
    nbelow=nασ_-Δασ_
    nabove=Δασ_
    # αbelow=solve_αX_from_nX(ebelow,nbelow,βbelow,nασ_)
    # αabove=solve_αX_from_nX(eabove,nabove,βabove,1-nασ_)
    nbelowcheck,Abelow,Kbelow,nkbelow=cal_nX_AX_KX_nkX_safe(ebelow,nbelow,βbelow,nασ_)
    nabovecheck,Aabove,Kabove,nkabove=cal_nX_AX_KX_nkX_safe(eabove,nabove,βabove,1-nασ_)
    Abelow,Aabove,Kbelow,Kabove,[nkbelow,nkabove]
end

"""
for compute Z
"""
function cal_αασ_(nασ_,Δασ_,βασ_,eασ_)
    βbelow,βabove=βασ_
    ebelow,eabove=eασ_
    nbelow=nασ_-Δασ_
    nabove=Δασ_
    αbelow=solve_αX_from_nX(ebelow,nbelow,βbelow,nασ_)
    αabove=solve_αX_from_nX(eabove,nabove,βabove,1-nασ_)
    [αbelow,αabove]
end
