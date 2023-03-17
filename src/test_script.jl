# band structure
# we move the data to es_files
e_fn=gene_spline_band("./es_files/es_inf.dat")
nσ=[0.4,0.6]
ϵsσ=[gene_ϵs(e_fn,nσ[1]),gene_ϵs(e_fn,nσ[2])]



