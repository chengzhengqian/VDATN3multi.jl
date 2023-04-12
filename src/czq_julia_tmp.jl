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


