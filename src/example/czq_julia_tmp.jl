for t1 in t1s
    results=[]
    para_dict=Dict()
    model_n3=create_model_two_band_t1_t2_half(t1,t2)
    set_para(model_n3,[0.37,0.37,0.6,0.6,0.1/t1,0.1,0.1])
    for U in Us
        print("processing t1=$(t1) U=$(U)\n")
        push!(results,solve_model(model_n3,U))
        para_dict[U]=get_para(model_n3)
    end
    saveData(results,"$(data_dir)/two_band_t1_$(t1)_t2_$(t2)_half_result.dat")
    save_object("$(data_dir)/two_band_t1_$(t1)_t2_$(t2)_half_para.jld2",para_dict)
end
