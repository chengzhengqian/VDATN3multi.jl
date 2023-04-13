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
