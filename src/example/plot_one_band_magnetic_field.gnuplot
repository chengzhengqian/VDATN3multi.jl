set terminal qt size 800,600  font ",20"


n3file="./data_one_band_half/result_n3.dat"
# data=[U,  Etotal,  Eloc, Ek,  nn[1], Δασ[1], Aασ_above[1], Aασ_below[1], αασ[1][1], βασ[1][1],  G12ασ[1], Zασ[1]]

n3fileB0="./data_one_band_magnetic_field/one_band_half_B_0_result_n3.dat"
    # data=[U,  Etotal,  Eloc, Ek,  nασ[1],nασ[2], nn[1], Δασ[1],Δασ[2], Aασ_above[1], Aασ_above[2], Aασ_below[1],Aασ_below[2], αασ[1][1],αασ[1][2],βασ[1][1],βασ[1][2],  G12ασ[1],G12ασ[2], Zασ[1], Zασ[2]]

# we first check the density
set key left bottom
set xrange [0:10]
set autoscale y
set xlabel "U/t"
set ylabel "{/Symbol D}n_{{/Symbol s}}"
plot n3fileB0 u ($1):($5-0.5) w l t  "spin up",  n3fileB0 u ($1):($6-0.5) w l t  "spin down"

set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_B_0_dn_U.png"
replot
unset terminal

set key left bottom
set xrange [0:10]
set autoscale y
set xlabel "U/t"
set ylabel "100*{/Symbol D}n  or  M"
set key right bottom
plot n3fileB0 u ($1):(($5/2+$6/2-0.5)*100) w l t  "100*{/Symbol D}n",  n3fileB0 u ($1):($5-$6) w l t  "M"

set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_B_0_dn_M_U.png"
replot
unset terminal

set xrange [0:10]
set yrange [0:0.25]
set xlabel "U/t"
set ylabel "d"
set key right top
plot n3fileB0 u ($1):($7)  t  "B=0", n3file  u 1:5 w l  t "n=0.5"

set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_B_0_fixed_check_d_U.png"
replot
unset terminal

n3fileB(U)="./data_one_band_magnetic_field/one_band_half_B_U_".U."_result_n3.dat"
    # [B,U,  Etotal,  Eloc, Ek,  nασ[1],nασ[2], nn[1], Δασ[1],Δασ[2], Aασ_above[1], Aασ_above[2], Aασ_below[1],Aασ_below[2], αασ[1][1],αασ[1][2],βασ[1][1],βασ[1][2],  G12ασ[1],G12ασ[2], Zασ[1], Zασ[2]]
dmftB(U)="./dmft_data/one_band_inf_mag_h_u_".U.".txt"

Us="1.0 2.0 3.0 4.0 8.0"
Usdmft="1.0 2.0 3.0 4.0"
set xrange [0:1.2]
set yrange [0:0.9]
set xlabel "B/t"
set ylabel "M"
set key right bottom
unset label
set label "U=1.0" font ",12" at 0.7,0.5
set label "U=2.0" font ",12" at 0.55,0.61
set label "U=3.0" font ",12" at 0.47,0.77
set label "U=4.0" font ",12" at 0.31,0.8
set label "U=8.0" font ",12" at 0.12,0.75

plot  NaN lc "black" t "DMFT(NRG)", NaN lc "red" t "N=3" ,  for [U in Usdmft] dmftB(U) u ($1):($2) w l lc "black" t "" ,for [U in Us] n3fileB(U) u ($1):($6-$7) w l lc "red" t ""


set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_half_inf_M_B.png"
replot
unset terminal

