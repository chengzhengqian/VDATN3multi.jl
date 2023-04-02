# first example, solve the problem at half-filling

include("../vdat.jl")

N_spin_orbital_=2
symmetry_=[[1,2]]
n_target_=[0.5]
e_fn=gene_spline_band("../es_files/es_inf.dat")
e_fns_=[e_fn]
U=1.0
interaction_=gene_interaction(U,0,N_spin_orbital_)
chemical_potential_=[]
model_n3=create_model(N_spin_orbital_,symmetry_,n_target_,
                      interaction_,chemical_potential_,e_fns_;
                      particle_hole_symmetric=true,N_time_step=3)

results=[]
Us=1.0:0.1:10.0
# U=1.0  model=model_n3   model_n3.obs

for U in Us
    print("processing U=$(U)\n")
    model_n3.interaction=gene_interaction(U,0,N_spin_orbital_)
    res=optimize(p->(set_para(model_n3,p);compute(model_n3)),get_para(model_n3))
    set_para(model_n3,res.minimizer)
    compute(model_n3;all_obs=true)
    @get_obs model_n3 Etotal Eloc Ek  nn Δασ Aασ_above Aασ_below αασ βασ  G12ασ Zασ
    data=[U,  Etotal,  Eloc, Ek,  nn[1], Δασ[1], Aασ_above[1], Aασ_below[1], αασ[1][1], βασ[1][1],  G12ασ[1], Zασ[1]]
    push!(results,data)
end

# save the result to disk
data_dir="./data_one_band_half"
mkpath(data_dir)
filename="result_n3.dat"
saveData(results, "$(data_dir)/$(filename)")



# Now, we can similarly run N=2, (i.e Gutzwiller Approximation)
model_n2=create_model(N_spin_orbital_,symmetry_,n_target_,
                      interaction_,chemical_potential_,e_fns_;
                      particle_hole_symmetric=true,N_time_step=2)


results=[]
Us=1.0:0.1:10.0

# U=2.0
# model_n2.interaction
for U in Us
    print("processing U=$(U)\n")
    model_n2.interaction=gene_interaction(U,0,N_spin_orbital_)
    res=optimize(p->(set_para(model_n2,p);compute(model_n2)),get_para(model_n2))
    set_para(model_n2,res.minimizer)
    compute(model_n2;all_obs=true)
    @get_obs model_n2 Etotal Eloc Ek  nn Zασ
    data=[U,  Etotal,  Eloc, Ek,  nn[1],Zασ[1]]
    push!(results,data)
end

filename="result_n2.dat"
saveData(results, "$(data_dir)/$(filename)")
