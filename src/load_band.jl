

"""
inf_fn=gene_spline_band("./es_inf.dat")
d2_fn=gene_spline_band("./es_2d.dat")
filename="./es_files/es_inf.dat"
"""
function gene_spline_band(filename;scale=1.0)
    data=reshape(loadData(filename),:)*scale
    N_sample=size(data)[1]
    index=collect(linspace(0,1,N_sample))
    return Spline1D(index,data)
end
# this is not useful for the landau level
# """
# p(ϵ<E)~E
# use to compute DOE
# p_e_fn=gene_spline_band_inv("$(phi_B_data)/phi_B_$(tag)_process.dat")
# filename="$(phi_B_data)/phi_B_$(tag)_process.dat"
# """
# function gene_spline_band_inv(filename)
#     data=reshape(loadData(filename),:)
#     N_sample=size(data)[1]
#     index=collect(linspace(0,1,N_sample))
#     return Spline1D(data,index)
# end

"""
e_fn=inf_fn
n=0.4
ϵs=gene_ϵs(e_fn,nσ[1])
nσ
ϵsσ=[gene_ϵs(e_fn,nσ[1]),gene_ϵs(e_fn,nσ[2])]
Here, we use integrate to improve the accuracy.
"""
function gene_ϵs(e_fn,n;N_samples=40,N_minimal=4)
    index_points_below=max(trunc(Int64,n*N_samples),N_minimal)
    index_points_above=max(trunc(Int64,(1-n)*N_samples),N_minimal)
    index_below=collect(linspace(0,n,index_points_below+1))
    energy_below=[(integrate(e_fn,index_below[i],index_below[i+1]))/(index_below[i+1]-index_below[i]) for i in 1:index_points_below]
    index_above=collect(linspace(n,1,index_points_above+1))
    energy_above=[(integrate(e_fn,index_above[i],index_above[i+1]))/(index_above[i+1]-index_above[i]) for i in 1:index_points_above]
    [energy_below,energy_above]
end


