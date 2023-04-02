# using CZQUtils
# from CZQUtils, so one does not need to download that package
linspace(start,stop,length)=range(start,stop=stop,length=length)

saveData(data,filename)= open(filename,"w") do io writedlm(io,data) end

loadData(filename)=readdlm(filename)


# from previous utils.jl
# previously, we use a dictionary to store the information for the purpose of automatic differentiation

"""
directly compute
Γ=1,..,2^(N_spin_orbital)
cal_Γασ(1,4) -> 0 0 0 0
cal_Γασ(2,4) -> 0 0 0 1
# this is consistent with direct product order
# previous, we use a dictionary to store the data, benchmark which one we should use
"""
function cal_Γασ(Γ,N_spin_orbital)
    reverse(digits(Γ-1,base=2,pad=N_spin_orbital))
end

# there are two ways to control the density, one is using the two particle fluctuation, the other is using the uniform transform

"""
using fluctuation to construct w2(w^2)
each η correspons a N(N>=2) density fluctuation that won't change the single particle density matrix
now, we need a conventon to generate all idx (combination with n>=2)
i.e for N_spin_orbital=4
we have 
[1,2] ...
[1,2,3] ...
[1,2,3,4]..
The order does not matter, as long as we are consistent through the calculation
For example
N_spin_orbital=3
ηToIdx=gene_idx(N_spin_orbital)
 [1, 2]
 [1, 3]
 [2, 3]
 [1, 2, 3]

"""
function gene_idx(N_spin_orbital)
    [idx for idx in combinations(1:N_spin_orbital) if length(idx)>=2]
end


"""
idx=[i,j,..k] the collection of sites the density correlation evolves
descripte the contribution of a single functuation (parametrized by idx)  to w2
For example
N_spin_orbital=2
idx=[1,2]
cal_Vη(N_spin_orbital,idx)
  0.25
 -0.25
 -0.25
  0.25
idx=[1]
cal_Vη(N_spin_orbital,idx)
 -0.5
 -0.5
  0.5
  0.5

"""
function cal_Vη(N_spin_orbital,idx)
    N_Γ=2^N_spin_orbital
    Vη=zeros(N_Γ)
    for Γ in 1:N_Γ
        Γασ=cal_Γασ(Γ,N_spin_orbital)
        Vη[Γ]=prod((Γασ.-0.5)[idx])
    end
    Vη
end

"""
η, the fluctuation index for 2-particle and beyond
Γ, atomic index
VΓη, N_Γ×N_η, map a vector in η space (i.e variational parameters )to Γ space (i.e w2)
ηToIdx, a array of N-particle functional index, i.e  [idx...], 
VΓη,ηToIdx=cal_VΓη_ηToIdx(2)
VΓη=>
  0.25
 -0.25
 -0.25
  0.25
ηToIdx =>
[ [1, 2]]


"""
function cal_VΓη_ηToIdx(N_spin_orbital)
    ηToIdx=gene_idx(N_spin_orbital)
    N_Γ=2^N_spin_orbital
    N_η=length(ηToIdx)
    VΓη=zeros(N_Γ,N_η)
    for i in 1:N_η
        VΓη[:,i]=cal_Vη(N_spin_orbital,ηToIdx[i])
    end
    VΓη,ηToIdx
end

# function to combine local representation to full representation

"""
extend Xmatασ[idx] to the full space, by assume idx in  spin orbital

"""
function cal_Xmatfull(pmatwασ,Xmatασ,idx::Integer)
    kron(pmatwασ[1:(idx-1)]...,Xmatασ[idx],pmatwασ[(idx+1):end]...)
end

"""
extend Xmatασ[idx1]*Xmatασ[idx2] to the full space, assuming idx1<idx2
"""
function cal_Xmatfull_(pmatwασ,Xmatασ,idx1::Integer,idx2::Integer)
    kron(pmatwασ[1:(idx1-1)]...,Xmatασ[idx1],pmatwασ[(idx1+1):(idx2-1)]...,Xmatασ[idx2],pmatwασ[(idx2+1):end]...)
end

"""
idx1, idx2 should be different, but idx1<idx2 is not required
supertype(Integer)
subtypes(Integer)
"""
function cal_Xmatfull(pmatwασ,Xmatασ,idx1::Integer,idx2::Integer)
    if(idx1<idx2)
        cal_Xmatfull_(pmatwασ,Xmatασ,idx1,idx2)
    elseif(idx1>idx2)
        cal_Xmatfull_(pmatwασ,Xmatασ,idx2,idx1)
    else
        error("idx1 $(idx1) idx2 $(idx2) should be different")
    end    
end

# there are some function to constraint x, we put them back later!!

"""
compute expetation values, w^2 is assumed to be normalizd,
pmatwασ is just array of identity matrix in w basis. But we leave room for future general code
!!! assuming the normalization
# if we working on the u rerensetation, pmat is not identity matrix
!! add the general function later
"""
function expt(w,Xmatfull)
    dot(w,Xmatfull*w)
end


"""
compute for g11, g12 etc
"""
function cal_Xασ(w,pmatwασ,Xmatwασ)
    N_spin_orbital=length(pmatwασ)
    [expt(w,cal_Xmatfull(pmatwασ,Xmatwασ,i)) for i in 1:N_spin_orbital]
end

# add symmetry
"""
using symmetry to speed up the calculation
w must has the given symmetry (not check)
symmetry=[[spinorb1,spinorb2,...]...]
For entry in symmetry, we assume we have the permutation symmetry within the group.


"""
function cal_Xασ(w,pmatwασ,Xmatwασ,symmetry)
    N_spin_orbital=length(pmatwασ)
    result=zeros(N_spin_orbital)
    for (idx,term) in enumerate(symmetry)
        result[term].=expt(w,cal_Xmatfull(pmatwασ,Xmatwασ,symmetry[idx][1]))
    end
    result                      
end


"""
cal self-energy, G11=G22=1/2, g11=g22=n, G12=-G21,g12=-g21
return [[s11,s12]...]
for example
cal_Slocασ([0.5,0.5],[0.4,0.4],[0.3,0.3])
 [0.9641319942611191, 0.2654232424677188]
 [0.9641319942611191, 0.2654232424677188]

"""
function cal_Slocασ(nασ,G12ασ,g12ασ)
    cal_sloc11sloc12.(nασ,G12ασ,g12ασ)
end

"""
r paramterize the self-energy

example
cal_rασ([[1.0,0.01]]) ≈ 0.0
cal_rασ([[0.01,1.0]]) ≈ 1.0
"""
function cal_rασ(Slocασ)
    [cal_rIns11s12(sloc11,sloc12)  for (sloc11,sloc12) in Slocασ]
end

"""
cασ=cal_cασ(nασ,Δασ,rασ)
compute c assuming G11=G22=1/2
used to update G12
"""
function cal_cασ(nασ,Δασ,rασ)
    cal_chalf.(nασ,Δασ,rασ)
end

"""
update G12 from momentum part
r,c parametrize the self-energy
"""
function cal_G12ασ(nασ,Δασ,rασ,cασ)
    cal_g012Incr.(nασ,Δασ,rασ,cασ)
end

"""
restrict the range of Δασ
n<=n-Δ
n>=Δ
 Δασ ∈ [0, min(n,1-n)]
charge transfer. !! we may need to use two different cutoff?
# this has problem in n->0 or n->1
restrict_Δασ(Δασ,nασ;cutoff=0)
"""
function restrict_Δασ(Δασ,nασ;cutoff=1e-4)
    clamp.(Δασ,cutoff,min.(nασ .- cutoff,1.0 .- nασ .- cutoff))
end

"""
restrict physical density
"""
function restrict_nασ(nασ;cutoff=1e-4)
    clamp.(nασ,cutoff,1.0 - cutoff)
end

"""
restrict self-energy
"""
function restrict_Slocασ(Slocασ;cutoff=1e-5)
    N_spin_orbital=length(Slocασ)
    Slocασ=[max.(Slocασ[i],cutoff) for i in 1:N_spin_orbital]
end

"""
compute charge fluctuation
remove the previous redundant nασ signature
"""
function cal_Δασ(g12ασ,Slocασ)
    N_spin_orbital=length(g12ασ)
    # [ cal_delta(g12ασ[idx],Slocασ[idx]...) for idx in 1:N_spin_orbital]
    # to fix the AD problem
    [ cal_delta(g12ασ[idx],Slocασ[idx][1],Slocασ[idx][2]) for idx in 1:N_spin_orbital]
end


# there are also some function to load data, see whether to include them later

# we put the part about symmetry here

"""
extend the para to full form
symmetry=[[group1_idx1...],[group2_idx1,...]]
# symmetry=[[1,2],[3,4]];
para=[value_for_group1,value_for_group2,...]
"""
function extend_with_symmetry(para,symmetry,N_spin_orbital)
    result=zeros(N_spin_orbital)
    for (idx,term) in enumerate(symmetry)
        result[term].=para[idx]
    end
    result
end

"""
extend G_para[G12_group1,...] to full form
we change the name from cal to extend
"""
function extend_G(G_para,symmetry,N_spin_orbital)
    extend_with_symmetry(G_para,symmetry,N_spin_orbital)
end

"""
extend β_para_below,β_para_above to full form
β_para_below=[3.0,4.0]
β_para_above=[3.1,4.1]
"""
function extend_β(β_para_below,β_para_above,symmetry,N_spin_orbital)
    β_para_below_ext=extend_with_symmetry(β_para_below,symmetry,N_spin_orbital)
    β_para_above_ext=extend_with_symmetry(β_para_above,symmetry,N_spin_orbital)
    [[β_para_below_ext[i],β_para_above_ext[i]] for i in 1:N_spin_orbital]
end

"""
generate momentum density point
update from previous version, so it takes e_fns now.
"""
function cal_eασ(e_fns,nασ,symmetry)
    N_spin_orbital=length(nασ)
    eασ=Vector{Any}(undef,N_spin_orbital)
    for (idx,term) in enumerate(symmetry)
        es=gene_ϵs(e_fns[idx],nασ[term[1]])
        for i in term
            eασ[i]=es
        end
    end
    eασ
end

"""
α is the orbital idx (start form 1) and σ is the spin indx (1,2) for spin up and own
"""
function get_idx(α,σ)
    (α-1)*2+σ
end


"""
return the largest integer that is smaller then N_average
"""
function absolute_floor(x)
    floor_x=floor(x)
    if(floor_x==x)
        x-1
    else
        floor_x
    end    
end

"""
return the smallest integer that is bigger then x
"""
function absolute_ceil(x)
    ceil_x=ceil(x)
    if(ceil_x==x)
        x+1
    else
        ceil_x
    end    
end


"""
compute Gfull (irreducible form)
"""
function cal_Gfull(nασ_,G12ασ_,g12ασ_,Aασ_below_,Aασ_above_,Slocασ_)
    gextraασ_=cal_glocReducedInA(Aασ_below_,Aασ_above_,Slocασ_[1],Slocασ_[2])
    Gfullασ_=cal_g0fullIngG12(nασ_,G12ασ_,g12ασ_,gextraασ_[1],gextraασ_[2],gextraασ_[3],gextraασ_[4])
end
