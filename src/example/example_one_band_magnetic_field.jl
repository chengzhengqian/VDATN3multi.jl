# fourth example, solve the one band model at given magnetic field

include("../vdat.jl")
# include("../compute.jl")
data_dir="./data_one_band_magnetic_field"
mkpath(data_dir)

function create_one_band_half_filling_magnetic_field(B,U,N_time_step_)
    N_spin_orbital_=2
    symmetry_=[[1],[2]]
    n_target_=[]
    e_fn=gene_spline_band("../es_files/es_inf.dat")
    e_fns_=[e_fn,e_fn]
    interaction_=gene_interaction(U,0,N_spin_orbital_)
    chemical_potential_=[(1,U/2+B),(2,U/2-B)]
    model=create_model(N_spin_orbital_,symmetry_,n_target_,
                       interaction_,chemical_potential_,e_fns_;
                       w_mode="free",cal_Eeff=cal_Eeff_U,N_Ueff=1,
                       N_time_step=N_time_step_)
end

function solve_n3_half_filling_magnetic_field(model_n3,U,B)
    model_n3.interaction=gene_interaction(U,0,model_n3.N_spin_orbital)
    model_n3.chemical_potential=[(1,U/2+B),(2,U/2-B)]
    res=optimize(p->(set_para(model_n3,p);compute(model_n3)),get_para(model_n3))
    set_para(model_n3,res.minimizer)
    compute(model_n3;all_obs=true)
    @get_obs model_n3 Etotal Eloc Ek nασ  nn Δασ Aασ_above Aασ_below  aασ bασ  G12ασ  Zασ
    data=[U,  Etotal,  Eloc, Ek,  nασ[1],nασ[2], nn[1], Δασ[1],Δασ[2], Aασ_above[1], Aασ_above[2], Aασ_below[1],Aασ_below[2], aασ[1][1],aασ[1][2],bασ[1][1],bασ[1][2],  G12ασ[1],G12ασ[2], Zασ[1], Zασ[2]]
end

# we first solve the problem with B=0
# U=1.0
Us=1.0:0.1:10.0
model_n3=create_one_band_half_filling_magnetic_field(0,Us[1],3)
results_half=[]
para_dict=Dict() # save the optimize variational parameters
set_para(model_n3,[0.37,0.37,0.7,0.7,0.7,0.7,0.33,0.33,0.67])
for U in Us
    print("processing $(U)\n")
    data=solve_n3_half_filling_magnetic_field(model_n3,U,0.0)
    push!(results_half,data)
    para_dict[U]=get_para(model_n3)
end


using JLD2
save_object("$(data_dir)/one_band_half_B_0_para.jld2",para_dict)
saveData(results_half,"$(data_dir)/one_band_half_B_0_result_n3.dat")
# para_dict=load_object("$(data_dir)/one_band_half_B_0_para.jld2")
#
U=1.0
UtoBs=Dict(1.0=>0.0:0.1:1.2, 2.0=>0.0:0.1:0.8,3.0=>0.0:0.05:0.5,4.0=>0.0:0.025:0.3, 8.0=>0.0:0.01:0.1)
Us=[1.0,2.0,3.0,4.0]
# Us=[8.0]
# get_para(model_n3)
for U in Us
    set_para(model_n3,para_dict[U])
    results=[]
    para_B_dict=Dict()
    for B in UtoBs[U]
        print("processing $(U) $(B)\n")
        data=solve_n3_half_filling_magnetic_field(model_n3,U,B)
        para_B_dict[B]=get_para(model_n3)
        push!(results,[B,data...])
    end
    save_object("$(data_dir)/one_band_half_B_U_$(U)_para.jld2",para_B_dict)
    saveData(results,"$(data_dir)/one_band_half_B_U_$(U)_result_n3.dat")
end




# debug for large magnetizaton, for check U=4.0
# U=4.0
# para_B_dict=load_object("$(data_dir)/one_band_half_B_U_$(U)_para.jld2")
# # start form B=0.25
# para=para_B_dict[0.25]
# set_para(model_n3,para)
# result=[]
# for B in 0.26:0.01:0.30
#     print("B is $(B)\n")
#     data=solve_n3_half_filling_magnetic_field(model_n3,U,B)
#     push!(results,[B,data...])
#     para_B_dict[B]=get_para(model_n3)
# end



# we found the para with problem as
# para_problem=[0.784082901099845,0.6877321375737462,0.3682046424655314,3.6550958632283534,2.825413265195028,-1.2379529554900333,12.758671520998323,-4.009087728387872,11.570054542605892]
# model=model_n3
# para=get_para(model_n3)


