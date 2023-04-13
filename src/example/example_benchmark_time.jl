# benchmark time
include("../vdat.jl")

N_spin_orbital_=2
N_spin_orbital_=4
N_spin_orbital_=6
N_spin_orbital_=8
N_spin_orbital_=10
N_spin_orbital_=12
N_spin_orbital_=14
N_spin_orbital_=16

symmetry_=[collect(1:N_spin_orbital_)]
n_target_=[0.5]
e_fn=gene_spline_band("../es_files/es_inf.dat")
e_fns_=[e_fn]
U=1.0
# interaction_=gene_interaction(U,0,N_spin_orbital_)
interaction_=[(1,2,U)]
chemical_potential_=[]
model_n3=create_model(N_spin_orbital_,symmetry_,n_target_,
                      interaction_,chemical_potential_,e_fns_;
                      particle_hole_symmetric=true,N_time_step=3
)
@time compute(model_n3)


N_spin_orbital_=2
N_spin_orbital_=4
N_spin_orbital_=6
N_spin_orbital_=8
N_spin_orbital_=10
N_spin_orbital_=12
N_spin_orbital_=14

symmetry_=[collect(1:N_spin_orbital_)]
n_target_=[0.5]
e_fn=gene_spline_band("../es_files/es_inf.dat")
e_fns_=[e_fn]
U=1.0
# interaction_=gene_interaction(U,0,N_spin_orbital_)
interaction_=[(1,2,U)]
chemical_potential_=[]
model_n2=create_model(N_spin_orbital_,symmetry_,n_target_,
                      interaction_,chemical_potential_,e_fns_;
                      particle_hole_symmetric=true,N_time_step=2
)
@time compute(model_n2)




