# sixth example, we add J to the previous example.

include("../vdat.jl")

"""
N_w_para_fixed=4
"""
function cal_w_fixed_two_band_half_with_J(neffασ,w_para)
    U1,U2,U3,U4=w_para
    Eeff=[cal_Eeff_fixed_two_band_half_with_J(cal_Γασ(i,4),U1,U2,U3,U4) for i in 1:2^4]
    exp.(-(Eeff.-minimum(Eeff)))
end

"""
U1,1,2
U2,3,4
U3,1,3;  2,4
U4, 1,4; 2,3;
"""
function cal_Eeff_fixed_two_band_half_with_J(Γασ,U1,U2,U3,U4)
    Δn1,Δn2,Δn3,Δn4=Γασ.-0.5
    U1*Δn1*Δn2+U2*Δn3*Δn4+U3*(Δn1*Δn3+Δn2*Δn4)+U4*(Δn1*Δn4+Δn2*Δn3)
end

function create_model_two_band_t1_t2_half(t1,t2)
    N_spin_orbital_=4
    symmetry_=[[1,2],[3,4]]
    n_target_=[0.5,0.5]
    e_fn_1=gene_spline_band("../es_files/es_inf.dat";scale=t1)
    e_fn_2=gene_spline_band("../es_files/es_inf.dat";scale=t2)
    e_fns_=[e_fn_1,e_fn_2]
    U=1.0
    interaction_=gene_interaction(U,0,N_spin_orbital_)
    chemical_potential_=[]
    model_n3=create_model(N_spin_orbital_,symmetry_,n_target_,
                          interaction_,chemical_potential_,e_fns_;
                          particle_hole_symmetric=true,N_time_step=3,
                          w_mode="fix",N_w_para_fixed=4,cal_w_fixed=cal_w_fixed_two_band_half_with_J)
end


"""
solve the model
"""
function solve_model(model_n3,U,J)
    model_n3.interaction=gene_interaction(U,J,N_spin_orbital_)
    res=optimize(p->(set_para(model_n3,p);compute(model_n3)),get_para(model_n3))
    set_para(model_n3,res.minimizer)
    compute(model_n3;all_obs=true)
    @get_obs model_n3 Etotal Eloc Ek  nn Δασ Aασ_above Aασ_below αασ βασ  G12ασ Zασ
    data=[U,  Etotal,  Eloc, Ek,  nn..., Δασ[[1,3]]..., Aασ_above[[1,3]]..., Aασ_below[[1,3]]..., G12ασ[[1,3]]..., Zασ[[1,3]]...]
end

using JLD2
N_spin_orbital_=4
data_dir="./two_band_half_t1_t2_with_J"
mkpath(data_dir)

t1=t2=1.0
Us=0.3:0.1:10.0
t1s=0.1:0.1:1.0
J_U=0.05
J_U=0.1


for t1 in t1s
    results=[]
    para_dict=Dict()
    model_n3=create_model_two_band_t1_t2_half(t1,t2)
    set_para(model_n3,[0.37,0.37,0.6,0.6,0.1/t1,0.1,0.1,0.1])
    for U in Us
        print("processing t1=$(t1) U=$(U)\n")
        push!(results,solve_model(model_n3,U,U*J_U))
        para_dict[U]=get_para(model_n3)
    end
    saveData(results,"$(data_dir)/two_band_t1_$(t1)_t2_$(t2)_J_U_$(J_U)_half_result.dat")
    save_object("$(data_dir)/two_band_t1_$(t1)_t2_$(t2)_J_U_$(J_U)_half_para.jld2",para_dict)
end



