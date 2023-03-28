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
