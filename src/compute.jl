# main part to create and compute model

# define options for Model (one could always define keys in model.options)
@define_option "require_precompute_eασ"
@define_option "is_density_fixed"
@define_option "w_mode"
@define_option "nασ_tolerance"
@define_option "cutoff_Sloc"
@define_option "cutoff_Δ"



"""
by default, we use fluctuation
add more option later
we first treat , w_mode=="fluc", add checking code later
# !! we update to fluc with symmetry
# we can also generalize to N_time_step 1-3
# add other mode later
add w_mode=="free"
# the problem is that we need to more flexible way
#  for w_mode=="free"
cal_Eeff=cal_Eeff_test,N_Ueff=0
# we now rename it from create_model_N3 to create_model
so we could include N=2
"""
function create_model(N_spin_orbital,symmetry,n_target,interaction,chemical_potential,e_fns;particle_hole_symmetric=false,w_mode="fluc",cal_Eeff=cal_Eeff_test,N_Ueff=0,N_time_step=3)
    N_symmetry=length(symmetry)
    # for place holder
    if(N_time_step==3)
        G12_para=[0.4 for i in 1:N_symmetry]
        if(particle_hole_symmetric)
            β_para=[1.0 for i in 1:N_symmetry]
        else
            β_para=[1.0 for i in 1:(2*N_symmetry)]
        end
    elseif(N_time_step==2)
        G12_para=[0.5 for i in 1:N_symmetry]
        β_para=[]               # for N=2, β=∞, we don't need
    else
        error("unsupported N_time_step $(N_time_step)")
    end

    options=Dict{String,Any}()
    obs=Dict{String,Any}()
    if(w_mode=="fluc")
        vΓηsym=VΓηsym(N_spin_orbital,symmetry)
        w_para=zeros(size(vΓηsym.VΓηsym)[2])
        options["vΓηsym"]=vΓηsym
        options["regulate_knorm"]=regulate_knorm_exp
        if(length(n_target)!=N_symmetry)
            error("w_mode is $(w_mode), so there should be $(N_symmetry) constraints on density, but n_target is $(n_target)!")
        end        
    elseif(w_mode=="free")
        options["N_Ueff"]=N_Ueff
        options["cal_Eeff"]=cal_Eeff
        w_para=zeros(N_symmetry+N_Ueff)
        if(length(n_target)>0)
            error("w_mode is $(w_mode), so there should be no constraint on density, but n_target is $(n_target)!")
        end        
    else
        error("w_mode $(w_mode) is not supported!")
    end
    model=Model(N_time_step,N_spin_orbital,symmetry,n_target,w_para,G12_para,β_para,interaction,chemical_potential,e_fns,options,obs)
    # set options
    # set density mode
    if(length(n_target)>0)
        # model.options
        option_is_density_fixed(model,true)
        option_require_precompute_eασ(model,true)
    else
        option_is_density_fixed(model,false)
        option_require_precompute_eασ(model,false)
    end
    option_w_mode(model,w_mode)
    # default options
    option_nασ_tolerance(model,1e-7)
    option_cutoff_Sloc(model,1e-5)
    option_cutoff_Δ(model,1e-5)
    model
end

# set and get variational parameters for the model

"""
first stage, set and get  variational parameters
para=get_para(model)
"""
function get_para(model::Model)
    if(model.N_time_step==3)
        [model.G12_para...,model.β_para...,model.w_para...]
    elseif(model.N_time_step==2)
        model.w_para
    else
        error("unsupported N_time_step $(N_time_step)")
    end    
end

"""
para=get_para(model)
set_para(model,para)
"""
function set_para(model::Model,para)
    if(model.N_time_step==3)
        N_G12=length(model.G12_para)
        N_β=length(model.β_para)
        N_w=length(model.w_para)
        model.G12_para[:]=para[1:N_G12]
        model.β_para[:]=para[(1+N_G12):(N_G12+N_β)]
        model.w_para[:]=para[(1+N_G12+N_β):(N_G12+N_β+N_w)]
    else
        model.w_para[:]=para
    end    
end

"""
comptue the discretization of k-points
"""
function compute_eασ(model::Model)
    model.obs["eασ"]=cal_eασ(model.e_fns,model.obs["nασ"],model.symmetry)
end


"""
compute the model
delete!(model.options,"nασ")
"""
function compute(model::Model)
    # now, we can start the calculaton, we first start with N=3, but add N=2 or N=1 later, model.options, model.obs, model.w_para
    treat_precompute_eασ(model)
    # for N=2,3, when N=2, compute_local already yield total energy
    compute_local(model)
    if(model.N_time_step==3)
        compute_momentum(model)
    end
    model.obs["Etotal"]
end


"""
for fixed density, we can precompute the discretization of k points
"""
function treat_precompute_eασ(model::Model)
    # if the density is fixed, and eασ has not been computed, we can precompute the discretization of the lattice, otherwise, we need to compute on the fly, or eασ is not valid.
    if(option_is_density_fixed(model))
        if(option_require_precompute_eασ(model))
            nασ=extend_with_symmetry(model.n_target,model.symmetry,model.N_spin_orbital)
            # model.obs["nασ"]=nασ
            @set_obs model nασ
            compute_eασ(model)
            # we don't need to recompute
            option_require_precompute_eασ(model,false)
        end
    end
end

"""
first stage, compute A block
"""
function compute_local(model::Model)
    # symmetry,N_spin_orbital=model.symmetry,model.N_spin_orbital
    # local symmetry,N_spin_orbital # force to check
    N_symmetry=length(model.symmetry)
    G12ασ=extend_G(model.G12_para,model.symmetry,model.N_spin_orbital)
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
        w=cal_w_free(w_para,model.symmetry,model.N_spin_orbital,cal_Eeff)
    else
        error("w_mode $(option_w_mode(model)) is not supported")
    end
    # compute g12, and check nασ
    g12ασ=cal_Xασ(w,pmatwασ,g12matwSασ,model.symmetry)
    nασ=cal_Xασ(w,pmatwασ,g11matwSασ,model.symmetry)
    # to fully check , set a option in model 
    # nασ=cal_Xασ(w,pmatwασ,g11matwSασ)
    # we need to check density constraint, or generate discretization of k points
    if(option_is_density_fixed(model))
        nασ_tolerance=option_nασ_tolerance(model)
        if(sum(abs.(model.obs["nασ"]-nασ))>nασ_tolerance)
            error("the density constraint error is larger than expected $(nασ_tolerance)")
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
        Δασ=restrict_Δασ(Δασ,nασ; cutoff=option_cutoff_Δ(model))
        # we put all assigment to obs here, model.obs
        @set_obs model G12ασ  w g12ασ Slocασ  Δασ pmatwασ g12matwSασ g11matwSασ
        # we can then compute the momentum part and BCD block
    elseif(model.N_time_step==2)
        # it is easier to just use GA
        # nασ=[0.2,0.3]
        g12ασ_0=sqrt.(nασ.*(1.0  .- nασ))
        Zασ=(g12ασ./g12ασ_0).^2
        nn=[expt(w,cal_Xmatfull(pmatwασ,g11matwSασ,idx1,idx2)) for (idx1,idx2,coefficient) in interaction]
        Eloc=sum([interaction[i][3]*nn[i]   for i in 1:length(interaction)])
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

"""
compute the momentum part and remaining blocks, and local interaction
model.chemical_potential=[(1,1.0)]
model.interaction
"""
function compute_momentum(model::Model)
    # model.obs
    @get_obs model nασ Δασ eασ G12ασ g12ασ Slocασ pmatwασ w
    N_symmetry=length(model.symmetry)
    if(length(model.β_para)==2*N_symmetry)
        β_below_para=model.β_para[1:N_symmetry]
        β_above_para=model.β_para[(1+N_symmetry):2*N_symmetry]
    elseif(length(model.β_para)==N_symmetry)
        β_below_para=β_above_para=model.β_para[1:N_symmetry]
    else
        error("length of model.β_para is $(length(model.β_para)), which should be either $(2*N_symmetry) or $(N_symmetry) ")
    end
    βασ=extend_β(β_below_para,β_above_para,model.symmetry,model.N_spin_orbital)
    Aασ_below,Aασ_above,αασ,nk,Ekασ=[Array{Any}(undef,model.N_spin_orbital) for _ in 1:5]
    # nασ,Δασ,βασ,eασ=model.obs["nασ"],model.obs["Δασ"],model.obs["βασ"],model.obs["eασ"]
    # nασ,Δασ,βασ,eασ=0,0,0,0, to test
    Ek=0 # kinetic energy
    # use the symmetry to simplify the calcuation
    for term in model.symmetry
        # just for debug, comment it later
        # global Ek,Aασ_below,Aασ_above,αασ,nk
        i=term[1]
        Aασ_below_,Aασ_above_,Kbelow_,Kabove_,αασ_,nk_=cal_Abelow_Aabove_Kbelow_Kabove_αασ_nk_(nασ[i],Δασ[i],βασ[i],eασ[i])
        Ekασ_=Kbelow_+Kabove_
        Ek+=Ekασ_*length(term)
        for idx in term
            Aασ_below[idx]=Aασ_below_
            Aασ_above[idx]=Aασ_above_
            αασ[idx]=αασ_
            nk[idx]=nk_
            Ekασ[idx]=Ekασ_
        end        
    end    
    # add to obs, model.obs, i=1, cal_Gfull -> 6 entries, G12, and remain 5 blocks
    g33matwασ=[cal_g33_mat_(cal_Gfull(nασ[i],G12ασ[i],g12ασ[i],Aασ_below[i],Aασ_above[i],Slocασ[i])) for i in 1:model.N_spin_orbital]
    nn=[expt(w,cal_Xmatfull(pmatwασ,g33matwασ,idx1,idx2)) for (idx1,idx2,coefficient) in model.interaction]
    Eloc=sum([model.interaction[i][3]*nn[i]   for i in 1:length(model.interaction)])
    # chemcial, notice model.chemical_potential maybe empty
    Echem=0.0
    for term in model.chemical_potential
        Echem+=-term[2]*nασ[term[1]]
    end
    Etotal=Ek+Eloc+Echem
    # compute Z , nασ=[0.5,0.5]
    # chemical potential for non-interacting density
    μασ=zeros(model.N_spin_orbital)
    for (idx,term) in enumerate(model.symmetry)
        i=term[1]
        μασ_=model.e_fns[idx](nασ[i])
        for ασ in term
            μασ[ασ]=μασ_
        end        
    end    
    Zασ=[cal_nk(αασ[i][1],βασ[i][1],μασ[i])-cal_nk(αασ[i][2],βασ[i][2],μασ[i]) for i in 1:model.N_spin_orbital]
    # model.obs
    @set_obs model  Etotal Ek Eloc Echem Aασ_below Aασ_above αασ  βασ nk nn g33matwασ Ekασ Zασ μασ 
    Etotal
end

