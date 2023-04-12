"""
neffασ=[0.5,0.5]
return the non-interacting w for neffασ
just a placeholder
"""
function cal_w_fixed_test(neffασ,w_para_fixed)
    w0=sqrt.(cal_w02(neffασ))
    return w0
end


# now, we need to develop a method to control the density

"""
w2 can be interprete as density, we fixed the density for spin-orbital idx to neff,
sum(w2)=1
w2=rand(4)
N_spin_orbital=2
idx=1
neff=0.4
neff=0.5
w2=rand(4)
w2=fix_density_for_w2(w2,1,0.8,N_spin_orbital)
w2=fix_density_for_w2(w2,2,0.4,N_spin_orbital)
w2=rand(2^4)
w2_check1=copy(w2)
w2_check1=fix_density_for_w2(w2_check1,1,0.8,N_spin_orbital)
w2_check1=fix_density_for_w2(w2_check1,2,0.4,N_spin_orbital)
w2_check2=copy(w2)
w2_check2=fix_density_for_w2(w2_check2,2,0.4,N_spin_orbital)
w2_check2=fix_density_for_w2(w2_check2,1,0.8,N_spin_orbital)
w2_check1-w2_check2
# it seems that the operation commutes with each other
# one can verify it using 4x4 matrix
"""
function fix_density_for_w2(w2,idx,neff,N_spin_orbital)
    w2=w2/sum(w2)
    idx_Γασp0,idx_Γασp1=gene_idx_Γασp0_idx_Γασp1(idx,N_spin_orbital)
    n=sum(w2[idx_Γασp1])
    if(n>neff)
        # transform from Γασp1 to Γασp0
        Δn=n-neff
        λ=Δn/n
        Δw2=w2[idx_Γασp1]*λ
        w2[idx_Γασp0].+=Δw2
        w2[idx_Γασp1].-=Δw2
        return w2
    elseif(n<neff)
        #transform from Γασp0 to Γασp1
        Δn=neff-n
        λ=Δn/(1-n)
        Δw2=w2[idx_Γασp0]*λ
        w2[idx_Γασp1].+=Δw2
        w2[idx_Γασp0].-=Δw2
        return w2
    else
        return w2
    end
end

"""
include("./vdat.jl")
for all the densities
w2=rand(4)
neffασ=[0.3,0.8]
w2=fix_density_w2(w2,neffασ)
"""
function fix_density_w2(w2,neffασ)
    N_spin_orbital=length(neffασ)
    for (idx,neff) in enumerate(neffασ)
        w2=fix_density_for_w2(w2,idx,neff,N_spin_orbital)
    end
    w2
end


"""
idx=6
Γασ=cal_Γασ(idx,4)
Γασ_to_idx(Γασ)-idx
"""
function Γασ_to_idx(Γασ)
    return sum(Γασ .* (2 .^ collect(length(Γασ)-1:-1:0)))+1
end

"""
idx=1
idx_Γασp0,idx_Γασp1=gene_idx_Γασp0_idx_Γασp1(idx,N_spin_orbital)
"""
function gene_idx_Γασp0_idx_Γασp1(idx,N_spin_orbital)
    N_spin_orbital_env=N_spin_orbital-1
    idx_Γασp0=[]
    idx_Γασp1=[]
    for i in 1:(2^N_spin_orbital_env)
        # state for environment   ;  i=1
        Γασp=cal_Γασ(i,N_spin_orbital_env)
        Γασp1=insert!(copy(Γασp),idx,1)
        push!(idx_Γασp1,Γασ_to_idx(Γασp1))
        Γασp0=insert!(copy(Γασp),idx,0)
        push!(idx_Γασp0,Γασ_to_idx(Γασp0))
    end
    idx_Γασp0,idx_Γασp1
end


