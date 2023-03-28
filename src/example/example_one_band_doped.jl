# second example, solve the problem at doped case,
# to compare with DMFT (NRG), we consider

include("../vdat.jl")

data_dir="./data_one_band_doped"
mkpath(data_dir)

ntotals=[0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98,0.99,0.995]


function create_one_band_doped(n_,U,N_time_step_)
    N_spin_orbital_=2
    symmetry_=[[1,2]]
    n_target_=[n_]
    e_fn=gene_spline_band("../es_files/es_inf.dat")
    e_fns_=[e_fn]
    interaction_=gene_interaction(U,0,N_spin_orbital_)
    chemical_potential_=[]
    create_model(N_spin_orbital_,symmetry_,n_target_,
                          interaction_,chemical_potential_,e_fns_;
                          particle_hole_symmetric=true,N_time_step=N_time_step_)
end


function solve_n3(model_n3,U)
    model_n3.interaction=gene_interaction(U,0,model.N_spin_orbital)
    res=optimize(p->(set_para(model_n3,p);compute(model_n3)),get_para(model_n3))
    set_para(model_n3,res.minimizer)
    compute(model_n3)
    @get_obs model_n3 Etotal Eloc Ek  nn Δασ Aασ_above Aασ_below αασ βασ  G12ασ Zασ
    data=[U,  Etotal,  Eloc, Ek,  nn[1], Δασ[1], Aασ_above[1], Aασ_below[1], αασ[1][1], βασ[1][1],  G12ασ[1], Zασ[1]]
end


Us=0.5:0.1:10.0
for ntotal_ in ntotals
    print("processing ntotal=$(ntotal_)\n")
    model_n3=create_one_band_doped(ntotal_/2,Us[1],3)
    set_para(model_n3,[0.4,1.0,-0.1])
    results=[]
    for U in Us
        print("processing U=$(U)\n")
        push!(results,solve_n3(model_n3,U))
    end
    saveData(results,"$(data_dir)/result_n3_ntotal_$(ntotal_).dat")
end    

# for n=2
function solve_n2(model_n2,U)
    model_n2.interaction=gene_interaction(U,0,model_n2.N_spin_orbital)
    res=optimize(p->(set_para(model_n2,p);compute(model_n2)),get_para(model_n2))
    set_para(model_n2,res.minimizer)
    compute(model_n2)
    @get_obs model_n2 Etotal Eloc Ek  nn Zασ
    data=[U,  Etotal,  Eloc, Ek,  nn[1],Zασ[1]]
end


for ntotal_ in ntotals
    print("processing ntotal=$(ntotal_)\n")
    model_n2=create_one_band_doped(ntotal_/2,Us[1],2)
    set_para(model_n2,[-0.1])
    results=[]
    for U in Us
        print("processing U=$(U)\n")
        push!(results,solve_n2(model_n2,U))
    end
    saveData(results,"$(data_dir)/result_n2_ntotal_$(ntotal_).dat")
end    
