# third example, solve the one band model at  given chemical potential
include("../vdat.jl")
# include("../compute.jl")
data_dir="./data_one_band_free_density"
mkpath(data_dir)



function create_one_band_free_density(Δμ,U,N_time_step_)
    N_spin_orbital_=2
    symmetry_=[[1,2]]
    n_target_=[]
    e_fn=gene_spline_band("../es_files/es_inf.dat")
    e_fns_=[e_fn]
    interaction_=gene_interaction(U,0,N_spin_orbital_)
    chemical_potential_=[(1,U/2+Δμ),(2,U/2+Δμ)]
    model=create_model(N_spin_orbital_,symmetry_,n_target_,
                       interaction_,chemical_potential_,e_fns_;
                       w_mode="free",cal_Eeff=cal_Eeff_U,N_Ueff=1,
                       N_time_step=N_time_step_)
end



function solve_n3_Δμ(model_n3,U,Δμ)
    model_n3.interaction=gene_interaction(U,0,model_n3.N_spin_orbital)
    model_n3.chemical_potential=[(1,U/2+Δμ),(2,U/2+Δμ)]
    res=optimize(p->(set_para(model_n3,p);compute(model_n3)),get_para(model_n3))
    set_para(model_n3,res.minimizer)
    compute(model_n3;all_obs=true)
    @get_obs model_n3 Etotal Eloc Ek nασ  nn Δασ Aασ_above Aασ_below  αασ βασ  G12ασ  Zασ
    # we temperarily remove Z and α
    data=[U,  Etotal,  Eloc, Ek,  nασ[1], nn[1], Δασ[1], Aασ_above[1], Aασ_below[1], αασ[1][1],βασ[1][1],  G12ασ[1], Zασ[1]]
end


# we first solve the problem with Δμ=0
Us=1.0:0.1:10.0
model_n3=create_one_band_free_density(0,Us[1],3)
set_para(model_n3,[0.37,1.0,1.0,0.5,1.0])
results_half=[]
para_dict=Dict() # save the optimize variational parameters
for U in Us
    print("processing $(U)\n")
    data=solve_n3_Δμ(model_n3,U,0.0)
    push!(results_half,data)
    para_dict[U]=get_para(model_n3)
end

using JLD2
save_object("$(data_dir)/one_band_dmu_0_para.jld2",para_dict)
saveData(results_half,"$(data_dir)/one_band_dmu_0_result_n3.dat")

# para_dict=load_object("$(data_dir)/one_band_dmu_0_para.jld2")
# now, we have checked the two working mode yield identical results at half-filling, we can use the free density model to probe how system responds the chemical potential.

Us=1.0:1.0:9.0
# we use the dmft result to estimate the range of Δμ we want to generate

UtoΔμs=Dict(1.0=>0.0:0.1:2.1,2.0=>0.0:0.1:2.6,
            3.0=>0.0:0.1:3.1,4.0=>0.0:0.1:3.5,
            5.0=>0.0:0.1:4.0,6.0=>0.0:0.1:4.5,
            7.0=>0.5:0.1:5.0,8.0=>1.0:0.1:5.5,
            9.0=>1.5:0.1:6.0,10.0=>2.0:0.1:6.5,
            )

# model=model_n3
# get_para(model); compute(model)
# Us=8.0:1.0:9.0
for U in Us
    para=para_dict[U]
    print("processing $(U) get para $(para)\n")
    set_para(model_n3,para)
    results=[]
    para_doped_dict=Dict()
    for Δμ in UtoΔμs[U]
        print("processing $(U) $(Δμ)\n")
        data=solve_n3_Δμ(model_n3,U,Δμ)
        para_doped_dict[Δμ]=get_para(model_n3)
        push!(results,[Δμ,data...])
    end
    save_object("$(data_dir)/one_band_dmu_U_$(U)_para.jld2",para_doped_dict)
    saveData(results,"$(data_dir)/one_band_U_$(U)_dmu_result_n3.dat")
end


# now, we run from some doped point to the half-filling


UtoΔμsReverse=Dict(9.0=>4.0:-0.1:1.5,8.0=>3.5:-0.1:1.0,7.0=>3.0:-0.1:0.0,6.0=>2.0:-0.1:0.0)
Us=[9.0,8.0,7.0,6.0]
Us=[6.0,7.0]
for U in Us
    para_doped_dict=load_object("$(data_dir)/one_band_dmu_U_$(U)_para.jld2")
    Δμs=UtoΔμsReverse[U]
    set_para(model_n3,para_doped_dict[Δμs[1]])
    results=[]
    for Δμ in Δμs
        print("processing $(U) $(Δμ)\n")
        data=solve_n3_Δμ(model_n3,U,Δμ)
        push!(results,[Δμ,data...])
    end
    saveData(results,"$(data_dir)/one_band_U_$(U)_dmu_result_n3_reverse.dat")
end
