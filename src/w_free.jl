# use exponetial (i.e effectic energy) to construct w,
# also, use the uniform transform the fixed density
# !! we allow a general density free way to construct interacting projector


# we first set up testing case
# we also test with symmetry
# N_spin_orbital=4
# symmetry=[[1,2],[3,4]]
# n_target=[0.5,0.6]
# G_para=[0.4,0.42]
# # we use ασ to indicates the full form, use target to indicate constraint, para indicate variaitional form
# nασ=extend_with_symmetry(n_target,symmetry,N_spin_orbital)
# G12ασ=extend_G(G_para,symmetry,N_spin_orbital)
# neffασ=cal_neffασ(nασ,G12ασ)
# p0=cal_w02(neffασ)
# build test distribution
# p=p0.*exp.(rand(length(p0)))
# we can use effective energy to build the general term, and respect the symmetry

# p=cal_p_proj(p0,Ueff_para,μeff_para,symmetry, N_spin_orbital)
# sum(p)
# we now need to first density, we do it per group of the symmetry

# """
# unpack Ueff_para to Ueff1,Ueff2,Ueff3
# unpack_Ueff([1.0])
# unpack_Ueff([1.0,0.2])
# unpack_Ueff([1.0,0.2,0.3])
# unpack_Ueff([1.0,0.2,0.3,4.0])

# """
# function  unpack_Ueff(Ueff_para)
#     n_Ueff=length(Ueff_para)
#     if(n_Ueff==1)
#         [Ueff_para[1],Ueff_para[1],Ueff_para[1]]
#     elseif(n_Ueff==2)
#         [Ueff_para[1],Ueff_para[2],Ueff_para[2]]
#     elseif(n_Ueff==3)
#         [Ueff_para[1],Ueff_para[2],Ueff_para[3]]
#     else
#         error("Ueff_para should have 1-3 elements")
#     end
# end
    
# """
# we present the general approach, 
# Ueff_para=[1.0]                 #  
# Ueff_para=[1.0,0.2]             # 
# Ueff_para=[1.0,0.2,0.1]         # 
# Ueff1 for the intraband interaction, 
# Ueff2 interband interaction with opposite spin, 
# Ueff3 and interband interaction with same spin

# μeff_para=[0.1,0.3]
# we can merge this to cal_w_xx in future
# we can use it either as free form or constraint density with the transform method
# """
# function cal_p_proj(p0,Ueff_para,μeff_para,symmetry, N_spin_orbital)
#     μeffασ=extend_with_symmetry(μeff_para,symmetry,N_spin_orbital)
#     N_orbital=trunc(Int,N_spin_orbital/2)
#     Ueff1,Ueff2,Ueff3=unpack_Ueff(Ueff_para)
#     # use three  effective interaction and chemical potential to parametrize the problem, p=p0*exp(-Eeff)
#     Eeff=[cal_Eeff_Ueff_μeff(Γ,N_spin_orbital,Ueff1,Ueff2,Ueff3,μeffασ) for Γ in 1:2^N_spin_orbital]
#     Eeff=Eeff.-minimum(Eeff)
#     p=p0.*exp.(-Eeff)
#     p=p/sum(p)
#     p
# end



# """
# the exponential parametrization, (we leave the density fixing part to future work, we will provide SU(N) half-filling case, and doped case)
# w_para=[μeff_para...,Ueff_para...]
# """
# function cal_w_exp(w_para,symmetry, N_spin_orbital)
#     N_symmetry=length(symmetry)
#     μeff_para=w_para[1:N_symmetry]
#     Ueff_para=w_para[(N_symmetry+1):end]
#     p0=cal_w02([0.5 for i in 1:N_spin_orbital ])
#     p=cal_p_proj(p0,Ueff_para,μeff_para,symmetry, N_spin_orbital)
#     sqrt.(p)
# end


"""
we compute w use an customize call_Eeff
the total w_para is N_symmetry (i.e. for μασ) + length(para)
cal_Eeff(Γασ,μασ,Ueff_para)
w_para=[0.1,0.2]
symmetry=[[1,2]]
cal_Eeff=cal_Eeff_test
w_para=[0.3,0.2,0.1,0.2]
symmetry=[[1,2]]
cal_Eeff=cal_Eeff_U_J
"""
function cal_w_free(w_para,symmetry,N_spin_orbital,cal_Eeff)
    N_orbital=trunc(Int,N_spin_orbital/2)
    N_symmetry=length(symmetry)
    μeff_para=w_para[1:N_symmetry]
    Ueff_para=w_para[(N_symmetry+1):end]
    μeffασ=extend_with_symmetry(μeff_para,symmetry,N_spin_orbital)
    Eeff=[cal_Eeff(cal_Γασ(Γ,N_spin_orbital),μeffασ,Ueff_para) for Γ in 1:2^N_spin_orbital]
    Eeff=Eeff.-minimum(Eeff)
    p=exp.(-Eeff)
    p=p/sum(p)
    sqrt.(p)
end



"""
simpliest way to construct the projector, using effective energy.
Γ=1...,2^N_spin_orbital
Γ=4
# assuming N_spin_orbital is even number
cal_Eeff_Ueff_μeff(Γ,N_spin_orbital,Ueff1,Ueff2,Ueff3,μeffασ)

# we need a better interface allow user to customize it
"""
function cal_Eeff_U_J(Γασ,μeffασ,Ueff_para)
    N_spin_orbital=length(μeffασ)
    N_orbital=trunc(Int,N_spin_orbital/2)
    Ueff1,Ueff2,Ueff3=Ueff_para
    Eeff=0
    # we first compute effective chemical potential, i=1
    for i in 1:N_spin_orbital
        ni=Γασ[i]
        Eeff+=-ni*μeffασ[i]
    end
    # compute local effective interaction, orb=1
    # intraband
    for orb in 1:N_orbital
        ni=Γασ[get_idx(orb,1)]
        nj=Γασ[get_idx(orb,2)]
        Eeff+=ni*nj*Ueff1
    end
    # interband, opposite spin
    for orb1 in 1:(N_orbital-1)
        for orb2 in (orb1+1):(N_orbital)
            ni=Γασ[get_idx(orb1,1)]
            nj=Γασ[get_idx(orb2,2)]
            Eeff+=ni*nj*Ueff2
            ni=Γασ[get_idx(orb1,2)]
            nj=Γασ[get_idx(orb2,1)]
            Eeff+=ni*nj*Ueff2
        end
    end
    # interband, same spin
    for orb1 in 1:(N_orbital-1)
        for orb2 in (orb1+1):(N_orbital)
            ni=Γασ[get_idx(orb1,1)]
            nj=Γασ[get_idx(orb2,1)]
            Eeff+=ni*nj*Ueff3
            ni=Γασ[get_idx(orb1,2)]
            nj=Γασ[get_idx(orb2,2)]
            Eeff+=ni*nj*Ueff3
        end
    end
    Eeff
end
"""
only U
cal_Eeff_U([1,1,0,0],[1,1,1,1],[0])
cal_Eeff_U([1,1,0,0],[0,0,0,0],[1])
"""
function cal_Eeff_U(Γασ,μeffασ,Ueff_para)
    N_spin_orbital=length(μeffασ)
    Ueff=Ueff_para[1]
    Eeff=0
    # we first compute effective chemical potential, i=1
    for i in 1:N_spin_orbital
        ni=Γασ[i]
        Eeff+=-ni*μeffασ[i]
    end
    for orb1 in 1:(N_spin_orbital-1)
        for orb2 in (orb1+1):N_spin_orbital
        ni=Γασ[orb1]
        nj=Γασ[orb2]
        Eeff+=ni*nj*Ueff
        end
    end    
    Eeff
end

# # now
# # p=p/sum(p)
# # we now fix the density for id=1 (symmetry idex)
# id=1
# region_idx=symmetry[id]
# region_n=nασ[region_idx[1]]
# region_size=length(region_idx)
# # use the previous notation
# N_max=region_size
# N_average=region_size*region_n
# # not necessary
# # !!!, # the idea should be that we first fix the environment state and uniformally adjust indiviudal one. Maybe we don't need to implement this general form.
# # and the code should be similar to the SUN case, (just not normalized)
# # maybe we should implemnt this laler
# # N_below=trunc(Int,absolute_floor(N_average))
# # N_above=trunc(Int,absolute_ceil(N_average))
# p_proj_below=0.0
# p_proj_above=0.0
# p_proj_equal=0.0
# N_below_expt=0.0                # we divide the normalization later
# N_above_expt=0.0 

# for Γ in 1:2^N_spin_orbital
#     # for initial test
#     global p_proj_below,p_proj_above,p_proj_equal,N_below_expt,N_above_expt
#     Γασ=cal_Γασ(Γ,N_spin_orbital)
#     N_region=sum(Γασ[region_idx])
#     if(N_region<N_average)
#         p_proj_below+=p[Γ]
#         N_below_expt+=p[Γ]*N_region
#     elseif (N_region>N_average)
#         p_proj_above+=p[Γ]
#         N_above_expt+=p[Γ]*N_region
#     else
#         p_proj_equal+=p[Γ]
#     end    
# end
# # apply normalization
# N_below_expt/=p_proj_below
# N_above_expt/=p_proj_above


# # and we uniformally transform the probability according to the type of Γ

# # test code
# pmatwασ,g12matwSασ=cal_p_g12_mat(G12ασ)
# g11matwSασ=cal_g11_mat(G12ασ)
# w=sqrt.(p)
# w0=sqrt.(p0)
# g11ασ=cal_Xασ(w,pmatwασ,g11matwSασ)
# g12ασ=cal_Xασ(w,pmatwασ,g12matwSασ)
# g11ασ_0=cal_Xασ(w0,pmatwασ,g11matwSασ)
# g12ασ_0=cal_Xασ(w0,pmatwασ,g12matwSασ)

"""
user should override this function
"""
function cal_Eeff_test(Γασ,μeffασ,Ueff_para)
    0.0
end

