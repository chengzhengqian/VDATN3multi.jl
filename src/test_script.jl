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


