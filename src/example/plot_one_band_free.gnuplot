set terminal qt size 800,600  font ",20"


n3file="./data_one_band_half/result_n3.dat"
# data=[U,  Etotal,  Eloc, Ek,  nn[1], Δασ[1], Aασ_above[1], Aασ_below[1], aασ[1][1], bασ[1][1],  G12ασ[1], Zασ[1]]
n3filedmu0="./data_one_band_free_density/one_band_dmu_0_result_n3.dat"
# data=[U,  Etotal,  Eloc, Ek,  nασ[1], nn[1], Δασ[1], Aασ_above[1], Aασ_below[1], aασ[1][1], bασ[1][1],  G12ασ[1], Zασ[1]]


# we first check the density
set xrange [0:10]
set autoscale y
set xlabel "U/t"
set ylabel "{/Symbol D}n"
plot n3filedmu0 u ($1):($5-0.5) w l t  "{/Symbol D}{/Symbol m}=0"


set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_dmu_0_dn_U.png"
replot
unset terminal


set xrange [0:10]
set yrange [0:0.25]
set xlabel "U/t"
set ylabel "d"
plot n3filedmu0 u ($1):($6)  t  "{/Symbol D}{/Symbol m}=0", n3file  u 1:5 w l  t "n=0.5"

set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_dmu_0_fixed_check_d_U.png"
replot
unset terminal

# now, we plot
n3filedmu(U)="./data_one_band_free_density/one_band_U_".U."_dmu_result_n3.dat"
dmftfile(U)="./dmft_data/one_band_inf_n_mu_u_".U.".txt"
# data=[ Δμ, U,  Etotal,  Eloc, Ek,  nασ[1], nn[1], Δασ[1], Aασ_above[1], Aασ_below[1], aασ[1][1], bασ[1][1],  G12ασ[1], Zασ[1]]


set xlabel "{/Symbol D}{/Symbol m}/t"
set ylabel "n"
set xrange [0:8]
set yrange [0.5:0.95]
set key right bottom
Us="1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0"
plot  NaN t "DMFT(NRG)" lc "black" , NaN t "N=3" lc "red",  for  [U in Us] dmftfile(U) u (-$1):(1.0-$2*0.5) t "" w l lc "black" , for  [U in Us ] n3filedmu(U) u ($1):($6)  w l    lc "red" t  ""

set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_n_dmu_from_half.png"
replot
unset terminal

# add runs from doped regime to half-filling case
Us_rev="9.0 8.0 7.0 6.0"
n3filedmurev(U)="./data_one_band_free_density/one_band_U_".U."_dmu_result_n3_reverse.dat"

plot  NaN t "DMFT(NRG)" lc "black" , NaN t "N=3" lc "pink",NaN t "N=3(reverse)" lc "red",  for  [U in Us] dmftfile(U) u (-$1):(1.0-$2*0.5) t "" w l lc "black" , for  [U in Us ] n3filedmu(U) u ($1):($6)  w l    lc "pink" t  "", for [U in Us_rev]   n3filedmurev(U) u ($1):($6)  w l    lc "red" t  ""

set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_n_dmu_from_half_reverse.png"
replot
unset terminal

