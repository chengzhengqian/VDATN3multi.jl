# main part to create and compute model

# define options for Model (one could always define keys in model.options)
@define_option "require_precompute_eασ"
@define_option "is_density_fixed"
@define_option "w_mode"
@define_option "nασ_tolerance"
@define_option "cutoff_Sloc"
@define_option "cutoff_Δ"
@define_option "N_k_samples"
@define_option "N_k_samples_min"



"""
by default, we use fluctuation
add more option later
we first treat , w_mode=="fluc", add checking code later
# !! we update to fluc with symmetry
# we can also generalize to N_time_step 1-3
# add other mode later
add w_mode=="free"
# the problem is that we need to have a more flexible way
#  for w_mode=="free"
cal_Eeff=cal_Eeff_test,N_Ueff=0
# we now rename it from create_model_N3 to create_model
so we could include N=2
# w_mode="fix", similar to "fluc", but uses customized function to construct w (notice in free, the user provides the function to compute u instead of w, as in u, there is no restriction on  effective density)
# we add options for controling the resolution of k points.
"""
function create_model(N_spin_orbital,symmetry,n_target,interaction,chemical_potential,e_fns;particle_hole_symmetric=false,w_mode="fluc",cal_Eeff=cal_Eeff_test,N_Ueff=0,N_time_step=3, cal_w_fixed=cal_w_fixed_test,N_w_para_fixed=0,N_k_samples=40,N_k_samples_min=4,nασ_tolerance=1e-7,cutoff_Sloc=1e-5,cutoff_Δ=1e-5)
    N_symmetry=length(symmetry)
    # for place holder
    if(N_time_step==3)
        G12_para=[0.4 for i in 1:N_symmetry]
        if(particle_hole_symmetric)
            b_para=[1.0 for i in 1:N_symmetry]
        else
            b_para=[1.0 for i in 1:(2*N_symmetry)]
        end
    elseif(N_time_step==2)
        G12_para=[0.5 for i in 1:N_symmetry]
        b_para=[]               # for N=2, β=∞, we don't need
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
    elseif(w_mode=="fix")
        w_para=zeros(N_w_para_fixed)
        options["N_w_para_fixed"]=N_w_para_fixed
        options["cal_w_fixed"]=cal_w_fixed
        if(length(n_target)!=N_symmetry)
            error("w_mode is $(w_mode), so there should be $(N_symmetry) constraints on density, but n_target is $(n_target)!")
        end        
    else
        error("w_mode $(w_mode) is not supported!")
    end
    model=Model(N_time_step,N_spin_orbital,symmetry,n_target,w_para,G12_para,b_para,interaction,chemical_potential,e_fns,options,obs)
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
    # check with model.options
    # nασ_tolerance=1e-7,cutoff_Sloc=1e-5,cutoff_Δ=1e-5
    option_nασ_tolerance(model,nασ_tolerance)
    option_cutoff_Sloc(model,cutoff_Sloc)
    option_cutoff_Δ(model,cutoff_Δ)
    option_N_k_samples(model,N_k_samples)
    option_N_k_samples_min(model,N_k_samples_min)
    model
end

# set and get variational parameters for the model

"""
first stage, set and get  variational parameters
para=get_para(model)
"""
function get_para(model::Model)
    if(model.N_time_step==3)
        [model.G12_para...,model.b_para...,model.w_para...]
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
        N_b=length(model.b_para)
        N_w=length(model.w_para)
        model.G12_para[:]=para[1:N_G12]
        model.b_para[:]=para[(1+N_G12):(N_G12+N_b)]
        model.w_para[:]=para[(1+N_G12+N_b):(N_G12+N_b+N_w)]
    else
        model.w_para[:]=para
    end    
end

"""
comptue the discretization of k-points
"""
function compute_eασ(model::Model)
    model.obs["eασ"]=cal_eασ(model.e_fns,model.obs["nασ"],model.symmetry;N_samples=option_N_k_samples(model),N_minimal=option_N_k_samples_min(model))
end


"""
compute the model
delete!(model.options,"nασ")
"""
function compute(model::Model; all_obs=false)
    # now, we can start the calculaton, we first start with N=3, but add N=2 or N=1 later, model.options, model.obs, model.w_para
    treat_precompute_eασ(model)
    # for N=2,3, when N=2, compute_local already yield total energy    # model.obs
    compute_local(model)
    if(model.N_time_step==3)
        compute_momentum(model)
    end
    if(all_obs)
        compute_Z(model)
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
get_para(model)
set_para(model,[ 0.6268414068801784, 0.19702119823595263, 0.6261874088728155, 3.8287379976382985, 0.7190212450683349])
# check the problem for n->0 or 1
model=model_n3
G12ασ=[0.1,0.4]
"""
function compute_local(model::Model)
    # symmetry,N_spin_orbital=model.symmetry,model.N_spin_orbital
    # local symmetry,N_spin_orbital # force to check
    N_symmetry=length(model.symmetry)
    G12ασ=extend_G(model.G12_para,model.symmetry,model.N_spin_orbital)
    # we add one more restriction, !! tweak the lower bound in future
    # !! set this as an option, 
    # G12ασ=clamp.(G12ασ,0.2,0.5)
    # this should be put to option, i.e how to restrict
    # !! only for N_time_step=3 !!, fix the bug for N=2
    if(model.N_time_step==3)
        G12ασ=restrict_G12ασ_.(G12ασ)
    end    
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

"""
compute the momentum part and remaining blocks, and local interaction
model.chemical_potential=[(1,1.0)]
model.interaction
"""
function compute_momentum(model::Model)
    # model.obs
    @get_obs model nασ Δασ eασ G12ασ g12ασ Slocασ pmatwασ w
    N_symmetry=length(model.symmetry)
    if(length(model.b_para)==2*N_symmetry)
        b_below_para=model.b_para[1:N_symmetry]
        b_above_para=model.b_para[(1+N_symmetry):2*N_symmetry]
    elseif(length(model.b_para)==N_symmetry)
        b_below_para=b_above_para=model.b_para[1:N_symmetry]
    else
        error("length of model.b_para is $(length(model.b_para)), which should be either $(2*N_symmetry) or $(N_symmetry) ")
    end
    bασ=extend_b(b_below_para,b_above_para,model.symmetry,model.N_spin_orbital)
    Aασ_below,Aασ_above,nk,Ekασ=[Array{Any}(undef,model.N_spin_orbital) for _ in 1:5]
    # nασ,Δασ,bασ,eασ=model.obs["nασ"],model.obs["Δασ"],model.obs["bασ"],model.obs["eασ"]
    # nασ,Δασ,bασ,eασ=0,0,0,0, to test
    Ek=0 # kinetic energy
    # use the symmetry to simplify the calcuation, term=[1,2], term=[1], term=[2]
    for term in model.symmetry
        # just for debug, comment it later
        # global Ek,Aασ_below,Aασ_above,aασ,nk
        i=term[1]
        Aασ_below_,Aασ_above_,Kbelow_,Kabove_,nk_=cal_Abelow_Aabove_Kbelow_Kabove_nk_safe(nασ[i],Δασ[i],bασ[i],eασ[i])
        Ekασ_=Kbelow_+Kabove_
        Ek+=Ekασ_*length(term)
        for idx in term
            Aασ_below[idx]=Aασ_below_
            Aασ_above[idx]=Aασ_above_
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
    # compute Z , nασ=[0.5,0.5] # we put in a differnet function, as we may need some special treat for α
    # chemical potential for non-interacting density
    # μασ=zeros(model.N_spin_orbital)
    # for (idx,term) in enumerate(model.symmetry)
    #     i=term[1]
    #     μασ_=model.e_fns[idx](nασ[i])
    #     for ασ in term
    #         μασ[ασ]=μασ_
    #     end        
    # end
    # we put this in a seperate function
    # Zασ=[cal_nk(aασ[i][1],bασ[i][1],μασ[i])-cal_nk(aασ[i][2],bασ[i][2],μασ[i]) for i in 1:model.N_spin_orbital]
    # model.obs
    @set_obs model  Etotal Ek Eloc Echem Aασ_below Aασ_above  bασ nk nn g33matwασ Ekασ 
    Etotal
end

"""
one can directly update interaction, chemical potential, but for fixed density mode, one just compute eασ once, so one need to set the flag for require_precompute_eασ. And eασ will be updated when compute(model) is called next time.
"""
function update_n_target(model::Model,n_target)
    if(option_is_density_fixed(model))
        model.n_target[:]=n_target
        option_require_precompute_eασ(model,true)
    else
        error("model should be in density fixed mode and then can be set with n_target!")
    end    
end

"""
for N=3, we move the part in 'compute' to calcualte Z here
model=model_n3
idx=1
"""
function compute_Z(model::Model)
    if(model.N_time_step==3)
        @get_obs model nασ Δασ bασ eασ
        μασ=zeros(model.N_spin_orbital)
        aασ=Array{Any}(undef,model.N_spin_orbital)
        Zασ=zeros(model.N_spin_orbital)
        for (idx,term) in enumerate(model.symmetry)
            i=term[1]
            μασ_=model.e_fns[idx](nασ[i])
            bασ_=bασ[i]
            eασ_=eασ[i]
            Δασ_=Δασ[i]
            nασ_=nασ[i]
            bbelow,babove=bασ_
            ebelow,eabove=eασ_
            nbelow=nασ_-Δασ_
            nabove=Δασ_
            abelow,nμbelow=cal_aX_nμX(ebelow,nbelow,bbelow,nασ_,μασ_)
            aabove,nμabove=cal_aX_nμX(eabove,nabove,babove,1-nασ_,μασ_)
            aασ_=[abelow,aabove]
            Zασ_=nμbelow-nμabove
            for ασ in term
                μασ[ασ]=μασ_
                aασ[ασ]=aασ_
                Zασ[ασ]=Zασ_
            end        
        end
        # one may need to add limiting cases. Do it later!
        # aασ=[cal_aασ_(nασ[i],Δασ[i],bασ[i],eασ[i]) for i in 1:model.N_spin_orbital]
        # Zασ=[cal_nk(aασ[i][1],bασ[i][1],μασ[i])-cal_nk(aασ[i][2],bασ[i][2],μασ[i]) for i in 1:model.N_spin_orbital]
        @set_obs model aασ Zασ
    end    
end
