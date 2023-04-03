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
