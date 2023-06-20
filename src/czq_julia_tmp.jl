function compute_local(model::Model)
    # symmetry,N_spin_orbital=model.symmetry,model.N_spin_orbital
    # local symmetry,N_spin_orbital # force to check
    N_symmetry=length(model.symmetry)
    G12ασ=extend_G(model.G12_para,model.symmetry,model.N_spin_orbital)
    # we add one more restriction, !! tweak the lower bound in future
    # !! set this as an option, 
    # G12ασ=clamp.(G12ασ,0.2,0.5)
    # this should be put to option, i.e how to restrict
    G12ασ=restrict_G12ασ_.(G12ασ)
    pmatwασ,g12matwSασ=cal_p_g12_mat(G12ασ)
    g11matwSασ=cal_g11_mat(G12ασ)
    # now, we need to use different scheme
    if(option_w_mode(model)=="fluc")
        # check w_fluc.jl
        # x should not be pure zero x=rand(11)*0.00001
        x=model.w_para 
        @get_obs model nασ
        vΓηsym=model.options["vΓηsym"]
        regulate_knorm=model.options["regulate_knorm"]
        w=cal_w_fluc(x,nασ,G12ασ,vΓηsym,regulate_knorm)
    elseif(option_w_mode(model)=="free")
        w_para=model.w_para
        cal_Eeff=model.options["cal_Eeff"]
        # we use Vwu to transform u to w (in u, there is no restricution and in w, the distribution should be bounded by neff which depends on G12
        u=cal_u_free(w_para,model.symmetry,model.N_spin_orbital,cal_Eeff)
        Vwu_mat=cal_Vwu_mat(G12ασ)
        Vwu=kron(Vwu_mat...)
        w=Vwu*u
        w=w/sqrt(sum(w.^2))
    elseif(option_w_mode(model)=="fix")
        w_para_fixed=model.w_para 
        @get_obs model nασ
        neffασ=cal_neffασ(nασ,G12ασ)
        cal_w_fixed=model.options["cal_w_fixed"]
        w=cal_w_fixed(neffασ,w_para_fixed)
        w=w/sqrt(sum(w.^2))
    else
        error("w_mode $(option_w_mode(model)) is not supported")
    end
    # compute g12, and check nασ
    g12ασ=cal_Xασ(w,pmatwασ,g12matwSασ,model.symmetry)
    nασ=cal_Xασ(w,pmatwασ,g11matwSασ,model.symmetry)
    # nασ=restrict_nασ(nασ)
    # for the free case, nασ may out of range and g12ασ may also have problems. one should compute the boundary of density, i.e, when G12 is not 1/2, w shoudl be within some range to yield physical values.
    # to fully check , set a option in model 
    # nασ=cal_Xασ(w,pmatwασ,g11matwSασ)
    # we need to check density constraint, or generate discretization of k points
    if(option_is_density_fixed(model))
        nασ_tolerance=option_nασ_tolerance(model)
        dens_error=sum(abs.(model.obs["nασ"]-nασ))
        if(dens_error>nασ_tolerance)
            error("the density constraint error is larger than expected $(nασ_tolerance) which is $(dens_error)")
        end        
    else
        @set_obs model nασ
        compute_eασ(model)
    end
    if(model.N_time_step==3)
        # compute self-energy
        Slocασ=cal_Slocασ(nασ,G12ασ,g12ασ)
        Slocασ=restrict_Slocασ(Slocασ,cutoff=option_cutoff_Sloc(model))
        # compute charge transfer
        Δασ=cal_Δασ(g12ασ,Slocασ)
        # we could not use cutoff here for Δασ, if it is close to 0, we use fix the prolbem in the momentum calculation; remove the option from constructor too!
        # Δασ=restrict_Δασ(Δασ,nασ; cutoff=option_cutoff_Δ(model))
        Δασ=restrict_Δασ(Δασ,nασ; cutoff=0.0)
        # we put all assigment to obs here, model.obs
        @set_obs model G12ασ  w g12ασ Slocασ  Δασ pmatwασ g12matwSασ g11matwSασ
        # we can then compute the momentum part and BCD block
    elseif(model.N_time_step==2)
        # it is easier to just use GA
        # nασ=[0.2,0.3]
        g12ασ_0=sqrt.(nασ.*(1.0  .- nασ))
        Zασ=(g12ασ./g12ασ_0).^2
        nn=[expt(w,cal_Xmatfull(pmatwασ,g11matwSασ,idx1,idx2)) for (idx1,idx2,coefficient) in model.interaction]
        Eloc=sum([model.interaction[i][3]*nn[i]   for i in 1:length(model.interaction)])
        Ek=0
        @get_obs model eασ
        for term in model.symmetry
            i=term[1]
            Ek+=mean(eασ[i][1])*nασ[i]*Zασ[i]*length(term)
        end    
        Echem=0.0
        for term in model.chemical_potential
            Echem+=-term[2]*nασ[term[1]]
        end
        Etotal=Ek+Eloc+Echem
        @set_obs model  Etotal Ek Eloc Echem nn w Zασ
    else
        error("unsupported N_time_step $(N_time_step)")
    end    
end
