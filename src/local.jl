"""
compute pmatwασ,g12matwSασ
G12ασ=[0.4,0.4]
pmatwασ,g12matwSασ=cal_p_g12_mat(G12ασ)
"""
function cal_p_g12_mat(G12ασ::Array)
    pmatwασ=cal_pmatw.(G12ασ)
    g12matwSασ=cal_g12matwS.(G12ασ)
    pmatwασ,g12matwSασ
end

"""
compute g11matwSασ
to double check nασ from model
"""
function cal_g11_mat(G12ασ::Array)
    cal_g11matwS.(G12ασ)    
end

"""
Gfull=[g012,g013,g023,g031,g032,g033]
"""
function cal_g33_mat_(Gfull)
    cal_g33matw(Gfull[1],Gfull[2],Gfull[3],Gfull[4],Gfull[5],Gfull[6])
end

"""
compute g33matw (not symmetric)
"""
function cal_g33_mat(Gfullασ)
    cal_g33_mat_.(Gfullασ)
end

"""
effective energy, use to compute w02 to have given nασ
"""
function cal_neffασ(nασ,G12ασ)
    cal_neffInn.(nασ,G12ασ)
end

"""
w02 entry for Γ with given density constraint, 
neffασ is computed from cal_neffασ(nασ,G12ασ)
peffασ=1-neffασ
"""
function cal_w02_Γ(Γ,N_spin_orbital,neffασ,peffασ)
    Γασ=cal_Γασ(Γ,N_spin_orbital)
    prod(Γασ.*neffασ+(1.0 .- Γασ).*peffασ)
end

"""
compute the non-interacting w02, with given (effective) denstiy constraint
cal_w02([0.5,0.5])
cal_w02([0.5,0.5,0.5,0.5])
neffασ is computed from cal_neffασ(nασ,G12ασ)
"""
function cal_w02(neffασ)
    peffασ=1.0  .- neffασ
    N_spin_orbital=length(neffασ)
    N_Γ=2^N_spin_orbital
    w02=[cal_w02_Γ(Γ,N_spin_orbital,neffασ,peffασ) for Γ in 1:N_Γ]
end


