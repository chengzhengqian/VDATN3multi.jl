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