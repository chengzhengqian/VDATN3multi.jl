# compute momentum part
## !! αX βX, to ensuere that the name of these Lagrange multiplier does not comflict with the orbital index α, and inverse temperature β, we rename them to a and b, but for the old part of the code, we won't change the names. so αX->aX, βX->bX
## we now change the naming from α,β -> a,b


"""
X represents either below or above, i.e, < or > in the paper
esX, the array to energy point in region X, assuming equal probability
aX, the energy where nkX=0.5
bX, controls the fluctuation, bX is bigger when the local interaction is stronger
weightX, the weight of region X, in, X as <, weightX is nασ, X as >, weightX is 1-nασ
"""
function cal_nX(esX,aX,bX,weightX)
    nkX=[ cal_nk(aX,bX,esX_) for esX_ in esX]
    nX=mean(nkX)*weightX
end

"""
find the aX
# n_mean should be within 0.0 and 1.0
# or nX∈[0,weightX]
# nX=constraint_nX(nX,weightX)
aX_guess is computed by assuming esX has only one effective data point
"""
function solve_aX_from_nX(esX,nX,bX,weightX)
    n_mean=nX/weightX
    e_mean=mean(esX)
    aX_guess=cal_alpha(e_mean,bX,n_mean)
    aX=find_zero(aX->cal_nX(esX,aX,bX,weightX)-nX,aX_guess,Order0())
end

"""
this compute all necessary quantities for a given region
nX,AX,KX,nkX
nX : density (has been weighted)
AX : charge fluctuation (has been weighted)
KX : kinetic energy (has been weighted)
nkX: momentum density distribution
"""
function cal_nX_AX_KX_nkX(esX,aX,bX,weightX)
    nkX=[ cal_nk(aX,bX,esX_) for esX_ in esX]
    nX=mean(nkX)*weightX
    AX=mean(sqrt.(nkX.*(1.0 .- nkX)))*weightX
    KX=mean(esX.*nkX)*weightX
    nX,AX,KX,nkX
end


"""
compute momentum information
Abelow,Aabove,Kbelow,Kabove,[abelow,aabove],[nkbelow,nkabove]

to break Aασ to Aασ_below and Aασ_above, so we can use the same form as the differential case. See ...v3.jl for details
nασ_,Δασ_,bασ_,eασ_=nασ[i],Δασ[i],bασ[i],eασ[i]
# debug when nασ->1
1-nασ_-Δασ_
"""
function cal_Abelow_Aabove_Kbelow_Kabove_aασ_nk_(nασ_,Δασ_,bασ_,eασ_)
    bbelow,babove=bασ_
    ebelow,eabove=eασ_
    nbelow=nασ_-Δασ_
    nabove=Δασ_
    abelow=solve_aX_from_nX(ebelow,nbelow,bbelow,nασ_)
    aabove=solve_aX_from_nX(eabove,nabove,babove,1-nασ_)
    nbelowcheck,Abelow,Kbelow,nkbelow=cal_nX_AX_KX_nkX(ebelow,abelow,bbelow,nασ_)
    nabovecheck,Aabove,Kabove,nkabove=cal_nX_AX_KX_nkX(eabove,aabove,babove,1-nασ_)
    Abelow,Aabove,Kbelow,Kabove,[abelow,aabove],[nkbelow,nkabove]
end

# there are some function to compute derivatives, check them later. see scheme1xx


# in the previous way, we first  assume the eασX has only one point and estimate the aX, then solve aX. This works well for normal value, but for n->0 or n->1, this has singularities. So in these cases, we could simply assume that nk is flat (or batter approximation?)
"""
we merge solving aX with cal_nX_AX_KX_nkX
# new version
esX=[1,2,3]
nX=1e-8
bX=4.0
weightX=0.1+1e-8
cal_nX_AX_KX_nkX_with_nX(esX,nX,bX,weightX)
"""
function cal_nX_AX_KX_nkX_with_nX(esX,nX,bX,weightX)
    aX=solve_aX_from_nX(esX,nX,bX,weightX)
    cal_nX_AX_KX_nkX(esX,aX,bX,weightX)
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
esX,nX,bX,weightX=eabove,nabove,babove,1-nασ_
!! we setup some cutoff to solve problem. When nX_mean ->0 or 1 within nX_mean_cutoff, we assume density is flat
"""
function cal_nX_AX_KX_nkX_safe(esX,nX,bX,weightX;nX_mean_cutoff=1e-6)
    if(weightX==0.0)
        return cal_nX_AX_KX_nkX_with_nX_mean(0.5,esX,0.0)
    else
        nX_mean=nX/weightX
        if(nX_mean<nX_mean_cutoff || nX_mean>1.0-nX_mean_cutoff)
            return cal_nX_AX_KX_nkX_with_nX_mean(nX_mean,esX,weightX)
        else
            return cal_nX_AX_KX_nkX_with_nX(esX,nX,bX,weightX)
        end        
    end    
end

"""
including the limiting case, we don't save aασ now, we can provide a function to compute it if the user wants.
nασ_,Δασ_,bασ_,eασ_=nασ[i],Δασ[i],bασ[i],eασ[i]
"""
function cal_Abelow_Aabove_Kbelow_Kabove_nk_safe(nασ_,Δασ_,bασ_,eασ_;nX_mean_cutoff=1e-6)
    bbelow,babove=bασ_
    ebelow,eabove=eασ_
    nbelow=nασ_-Δασ_
    nabove=Δασ_
    # abelow=solve_aX_from_nX(ebelow,nbelow,bbelow,nασ_)
    # aabove=solve_aX_from_nX(eabove,nabove,babove,1-nασ_)
    nbelowcheck,Abelow,Kbelow,nkbelow=cal_nX_AX_KX_nkX_safe(ebelow,nbelow,bbelow,nασ_;nX_mean_cutoff=nX_mean_cutoff)
    nabovecheck,Aabove,Kabove,nkabove=cal_nX_AX_KX_nkX_safe(eabove,nabove,babove,1-nασ_;nX_mean_cutoff=nX_mean_cutoff)
    Abelow,Aabove,Kbelow,Kabove,[nkbelow,nkabove]
end

"""
for compute Z
"""
function cal_aασ_(nασ_,Δασ_,bασ_,eασ_)
    bbelow,babove=bασ_
    ebelow,eabove=eασ_
    nbelow=nασ_-Δασ_
    nabove=Δασ_
    abelow=solve_aX_from_nX(ebelow,nbelow,bbelow,nασ_)
    aabove=solve_aX_from_nX(eabove,nabove,babove,1-nασ_)
    [abelow,aabove]
end


"""
as in some cases, nX->0 or 1, 
esX,nX,bX,weightX=ebelow,nbelow,bbelow,nασ_
esX,nX,bX,weightX=eabove,nabove,babove,1-nασ_
and in this case, we may just return aX 0.0 (maybe Inf to simplfy the problem saving to txt)
μX=0.0
"""
function cal_aX_nμX(esX,nX,bX,weightX,μX;nX_mean_cutoff=1e-6)
    if(weightX==0.0)
        return (0.0, 0.5)
    else
        nX_mean=nX/weightX
        if(nX_mean<nX_mean_cutoff || nX_mean>1.0-nX_mean_cutoff)
            return (0.0, nX_mean)
        else
            aX=solve_aX_from_nX(esX,nX,bX,weightX)
            nμX=cal_nk(aX,bX,μX)
            return (aX,nμX)
        end
    end
end
