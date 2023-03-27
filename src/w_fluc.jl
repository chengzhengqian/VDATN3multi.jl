# using fluctuation to parametrization
# we should consider the symmetry, also, allow w_para to be zero
# we add symmetry
struct VΓηsym
    # matrix from ηsym to Γ (w2)
    VΓηsym
    # physical meaning of ηsym
    ηsymToIdx
    # N_spin_orbital
end

function ηIdxToηsym(ηToIdx_,symmetry,symmetry_dict)
    results=[0 for _ in 1:length(symmetry)]
    for idx in ηToIdx_
        results[symmetry_dict[idx]]+=1
    end
    results
end

"""
xtow2=Vxtow2(2)
we need to consider symmetry
ηToIdx_=ηToIdx[2]
# vΓηsym=VΓηsym(model::Model)    
"""
function VΓηsym(N_spin_orbital::Integer,symmetry)    
    VΓη,ηToIdx=cal_VΓη_ηToIdx(N_spin_orbital)
    symmetry_dict=Dict()
    for (idx,term) in enumerate(symmetry)
        # global symmetry_dict
        for term_ in term
            symmetry_dict[term_]=idx
        end        
    end
    # we just need to store the vector
    ηsym_dict=Dict()
    ηsymToIdx=[]
    # ηToIdx_= ηToIdx[1]  idx_=1
    for (idx_,ηToIdx_) in enumerate(ηToIdx)
        ηsym=ηIdxToηsym(ηToIdx_,symmetry,symmetry_dict)
        if(haskey(ηsym_dict,ηsym))
            ηsym_dict[ηsym]+=VΓη[:,idx_]
        else
            ηsym_dict[ηsym]=VΓη[:,idx_]
            push!(ηsymToIdx,ηsym)
        end        
    end
    # should we normalize the vector?
    VΓηsym_=zeros(2^(N_spin_orbital),length(ηsymToIdx))
    for (idx_, ηsymToIdx_) in enumerate(ηsymToIdx)
        # VΓηsym_[:,idx_]=ηsym_dict[ηsymToIdx_]
        VΓηsym_[:,idx_]=cal_k0_knorm(ηsym_dict[ηsymToIdx_])[1]
    end
    VΓηsym(VΓηsym_,ηsymToIdx)
end

# function Vxtow2(N_spin_orbital::Integer)
#     VΓη,ηToIdx=cal_VΓη_ηToIdx(N_spin_orbital)
#     Vxtow2(VΓη,ηToIdx,N_spin_orbital)
# end

function norm1(k)
    sum(abs.(k))
end


"""
k0 is the unit vector, knorm is the norm of k
"""
function cal_k0_knorm(k)
    knorm=norm1(k)
    k0=k/knorm
    k0,knorm
end


"""
used to constraint k
# it is acutally the scale, and with the new scheme, this is the only function we need
"""
function cal_maxnorm(w02,k0)
    kmax=minimum([-w02[i]/k0[i] for i in 1:length(w02)  if k0[i]<0])
end


"""
xtow2, store the information about how x transform to w2
regulate_knorm, take knorm,kmax and returns a scale smaller than kmax
knormscaled=regulate_knorm(knorm,kmax)
k=VΓη*x, in w2 space
example to test this function:
N_spin_orbital=4
nασ=[0.2,0.3,0.4,0.5]
G12ασ=[0.4,0.5,0.3,0.45]
xtow2=Vxtow2(N_spin_orbital)
pmatwασ,g12matwSασ=cal_p_g12_mat(G12ασ)
g11matwSασ=cal_g11_mat(G12ασ)
x=rand(2^N_spin_orbital-N_spin_orbital-1)*0.01
w=cal_w_fluc(x,nασ,G12ασ,xtow2,regulate_knorm_linear)
g11ασ=cal_Xασ(w,pmatwασ,g11matwSασ)
g12ασ=cal_Xασ(w,pmatwασ,g12matwSασ)
#equivalent to 
g11ασ=[expt(w,cal_Xmatfull(pmatwασ,g11matwSασ,i)) for i in 1:N_spin_orbital]
g12ασ=[expt(w,cal_Xmatfull(pmatwασ,g12matwSασ,i)) for i in 1:N_spin_orbital]
this is useful to constraint for small system
x could be view as w para
w=sqrt.(w02)# 
"""
# function cal_w_fluc(x,nασ,G12ασ,xtow2::Vxtow2,regulate_knorm)
#     neffασ=cal_neffασ(nασ,G12ασ)
#     w02=cal_w02(neffασ)
#     k=xtow2.VΓη*x
#     k0,knorm=cal_k0_knorm(k)
#     kmax=cal_maxnorm(w02,k0)
#     # we regulate the norm
#     knormscaled=regulate_knorm(knorm,kmax)
#     k_reg=k0*knormscaled
#     w=sqrt.(k_reg+w02)
# end

"""
the way treat symmetry
regulate_knorm=regulate_knorm_exp
xsym=rand(6)*0.01
nασ=[0.1,0.1,0.6,0.6]
G12ασ=[0.4,0.4,0.3,0.3]
pmatwασ,g12matwSασ=cal_p_g12_mat(G12ασ)
g11matwSασ=cal_g11_mat(G12ασ)
w=cal_w_func(xsym,nασ,G12ασ,vΓηsym,regulate_knorm_exp)
w=cal_w_func(zeros(6),nασ,G12ασ,vΓηsym,regulate_knorm_exp)
g11ασ=cal_Xασ(w,pmatwασ,g11matwSασ)
g12ασ=cal_Xασ(w,pmatwασ,g12matwSασ)

"""
function cal_w_fluc(xsym,nασ,G12ασ,vΓηsym::VΓηsym,regulate_knorm)
    neffασ=cal_neffασ(nασ,G12ασ)
    w02=cal_w02(neffασ)
    if(norm1(xsym)==0.0)
        return sqrt.(w02)
    else
        k=vΓηsym.VΓηsym*xsym
        k0,knorm=cal_k0_knorm(k)
        kmax=cal_maxnorm(w02,k0)
        knormscaled=regulate_knorm(knorm,kmax)
        k_reg=k0*knormscaled
        return sqrt.(k_reg+w02)
    end    
end


# there are various choice of choose regulate_form
"""
regulate_knorm_linear(0.0,0.2), exponetial decay back to zero after knorm>kmax
"""
function regulate_knorm_linear(knorm,kmax)
    factor=(1-1e-4)
    if(knorm>kmax*factor)
        kmax*factor*exp(-(knorm-kmax*factor)/kmax)
    else
        knorm
    end    
end

"""
for the easy of benchmarking, direct cutoff
"""
function regulate_knorm_linear_flat(knorm,kmax)
    factor=(1-1e-4)
    if(knorm>kmax*factor)
        kmax*factor
    else
        knorm
    end    
end


"""
expoential
regulate_knorm_exp(100,10)
"""
function regulate_knorm_exp(knorm,kmax)
    kmax*(1-exp(-knorm/kmax))
end


# we can add other choice of regulation function, do it later.

