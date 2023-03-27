# the starting point, solving the one band Hubbard model


include("../vdat.jl")
# include("../compute.jl")

using Optim
# half-filling

# load the band structure
N_spin_orbital_=2
symmetry_=[[1,2]]
n_target_=[0.5]
e_fn=gene_spline_band("../es_files/es_inf.dat")
e_fns_=[e_fn]
U=1.0
interaction=gene_interaction(U,0,N_spin_orbital_)
chemical_potential=[]
model=create_model(N_spin_orbital_,symmetry_,n_target_,
                      interaction,chemical_potential,e_fns_;
                      particle_hole_symmetric=true,N_time_step=3)
compute(model)




res=optimize(p->(set_para(model,p);compute(model)),get_para(model))
set_para(model,res.minimizer)
compute(model)
model.obs
get_para(model)
model.obs["nn"]


model2=create_model(N_spin_orbital_,symmetry_,n_target_,
                      interaction,chemical_potential,e_fns_;
                      particle_hole_symmetric=true,N_time_step=2)
compute(model2)
model2.obs
res=optimize(p->(set_para(model2,p);compute(model2)),get_para(model2))
set_para(model2,res.minimizer)
compute(model2)
model2.obs
set_para(model2,[-0.1])
# it seems that N=2 is fine


# fixed density multiband
N_spin_orbital_=4
symmetry_=[[1,2],[3,4]]
n_target_=[0.2,0.6]
e_fn=gene_spline_band("../es_files/es_inf.dat")
e_fns_=[e_fn,e_fn]
U=1.0
interaction=gene_interaction(U,0,N_spin_orbital_)
chemical_potential=[]
model=create_model_N3(N_spin_orbital_,symmetry_,n_target_,
                      interaction,chemical_potential,e_fns_;
                      particle_hole_symmetric=false)
compute(model)

res=optimize(p->(set_para(model,p);compute(model)),get_para(model))
set_para(model,res.minimizer)
compute(model)
model.obs
get_para(model)
v=model.options["vΓηsym"];
size(v.VΓηsym)
v.ηsymToIdx




# U=10.0
# model.interaction=gene_interaction(U,0,N_spin_orbital)
# free density case
e_fn=gene_spline_band("../es_files/es_inf.dat")
N_spin_orbital_=2
symmetry_=[[1,2]]
e_fns_=[e_fn]
n_target_=[]
U=1.0
interaction=gene_interaction(U,0,N_spin_orbital_)
Δμ=0.1
chemical_potential=[(1,U/2+Δμ),(2,U/2+Δμ)]
# chemical_potential=[(1,U)]
model=create_model_N3(N_spin_orbital_,symmetry_,n_target_,
                      interaction,chemical_potential,e_fns_;
                      w_mode="free",cal_Eeff=cal_Eeff_U,N_Ueff=1)

compute(model)
model.obs
get_para(model)
set_para(model,[0.4,1.0,1.0,0.1,0.2])

res=optimize(p->(set_para(model,p);compute(model)),get_para(model))
para=res.minimizer
set_para(model,para)
compute(model)
model.obs
