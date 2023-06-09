# band structure
# we move the data to es_files
e_fn=gene_spline_band("./es_files/es_inf.dat")
nσ=[0.4,0.6]
ϵsσ=[gene_ϵs(e_fn,nσ[1]),gene_ϵs(e_fn,nσ[2])]

# test model

e_fn=gene_spline_band("./es_files/es_inf.dat")
N_spin_orbital=4
symmetry=[[1,2],[3,4]]
n_target=[0.3,0.4]
e_fns=[e_fn,e_fn]
interaction=gene_interaction(U,J,N_spin_orbital)
chemical_potential=[]
model=create_model_N3(N_spin_orbital,symmetry,n_target,interaction,chemical_potential,e_fns)
model.options
model.obs
compute(model)
set_para(model,[0.34,0.41,1.1,1.2,1.3,1.4,[0.01 for i in 1:11]...])
get_para(model)
# we put the code to model.jl and compute.jl
model.obs
compute_momentum(model)

# !!!! after model.jl
# w_fluc
include("./vdat.jl")
# include("./compute.jl")
N_spin_orbital=4
symmetry=[[1,2],[3,4]]
n_target=[0.3,0.4]
e_fn=gene_spline_band("./es_files/es_inf.dat")
e_fns=[e_fn,e_fn]
U,J=1.0,0.1
interaction=gene_interaction(U,J,N_spin_orbital)
chemical_potential=[]
model=create_model_N3(N_spin_orbital,symmetry,n_target,interaction,chemical_potential,e_fns)
model.options
compute(model)
model.obs["Ekασ"]
# we update w_fluc
# include("./compute.jl")
# now, we add w_free
N_spin_orbital=4
symmetry=[[1,2],[3,4]]
n_target=[]
e_fn=gene_spline_band("./es_files/es_inf.dat")
e_fns=[e_fn,e_fn]
U,J=1.0,0.1
interaction=gene_interaction(U,J,N_spin_orbital)
chemical_potential=[(1,1.0),(3,1.0)]
model=create_model_N3(N_spin_orbital,symmetry,n_target,interaction,chemical_potential,e_fns;w_mode="free",cal_Eeff=cal_Eeff_U_J,N_Ueff=3)
model.options
model.obs
model.options
compute(model)

# consider the N=2


# set_precompute_eασ(model,true)

# we use macro the generate here
# """
# check whether we need to compute momemtum points
# # we should also initialize this option
# """
# function is_eασ_computed(model::Model)
#     if(haskey(model.options, "k_points_precomputed"))
#         model.options["k_points_precomputed"]
#     else
#         false
#     end    
# end

# """
# set whether we need to compute momemtum points
# """
# function set_is_eασ_computed(model::Model,val::Bool)
#     model.options["k_points_precomputed"]=val
# end

# function is_density_fixed(model::Model)
#     model.options["is_density_fixed"]=true
# end

# function get_is_density_fixed(model::Model)
#     model.options["is_density_fixed"]
# end


# use

# check fix density mode, using

include("./vdat.jl")

function cal_w_fixed_one_band_half(neffασ,w_para)
    Ueff=w_para[1]
    [exp(-Ueff),1,1,exp(-Ueff)]
end


N_spin_orbital_=2
symmetry_=[[1,2]]
n_target_=[0.5]
e_fn=gene_spline_band("./es_files/es_inf.dat")
e_fns_=[e_fn]
U=1.0
interaction_=gene_interaction(U,0,N_spin_orbital_)
chemical_potential_=[]
model_n3=create_model(N_spin_orbital_,symmetry_,n_target_,
                      interaction_,chemical_potential_,e_fns_;
                      particle_hole_symmetric=true,N_time_step=3,w_mode="fix",N_w_para_fixed=1,cal_w_fixed=cal_w_fixed_one_band_half)

model_n3.options
model=model_n3
model.obs
model.options
get_para(model_n3)
set_para(model_n3,[0.37,1.0,0.5])
compute(model_n3)
model_n3.obs

# 
# test for new options


# include("../src/vdat.jl")
include("./vdat.jl")
# check for "/ssh:gin:/burg/ccce/users/zc2255/VDATN3multi.jl/gin_run/SU_N_doped_n3_taylor.jl"

N_spin_orbital_=5
symmetry_=[collect(1:N_spin_orbital_)]
n=0.5
n_target_=[n]
e_fn_=gene_spline_band("../src/es_files/es_inf.dat";scale=1.0)
e_fns_=[e_fn_]
U=1.0
# one should use 
# interaction_=gene_interaction(U,0,N_spin_orbital_)
interaction=gene_interaction_degenerate_J_0(U,N_spin_orbital_)
chemical_potential_=[]

function cal_Eeff_SU_N_taylor(Γασ,paras)
    N_particle=sum(Γασ)
    μeff=paras[1]
    coeff=paras[2:end]
    sum([c*(N_particle-μeff)^i for (i,c) in enumerate(coeff)])
end

function cal_w_SU_N_taylor(neffασ,w_para)
    Eeff=[cal_Eeff_SU_N_taylor(cal_Γασ(i,N_spin_orbital_),w_para) for i in 1:2^N_spin_orbital_]
    w2=exp.(-(Eeff.-minimum(Eeff)))
    w2=fix_density_w2(w2,neffασ)
    sqrt.(w2/sum(w2))
end

# add more options
model_n3=create_model(N_spin_orbital_,symmetry_,n_target_,
                      interaction_,chemical_potential_,e_fns_;
                      particle_hole_symmetric=false,N_time_step=3,
                      w_mode="fix",N_w_para_fixed=N_spin_orbital_-1,cal_w_fixed=cal_w_SU_N_taylor,N_k_samples=100,N_k_samples_min=20,nασ_tolerance=1e-8,cutoff_Sloc=1e-6,cutoff_Δ=1e-6)

model_n3.options
model=model_n3
compute(model)
model.obs["eασ"][1][1]






