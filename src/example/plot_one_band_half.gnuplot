set terminal qt size 800,600  font ",20"

n3file="./data_one_band_half/result_n3.dat"
# data=[U,  Etotal,  Eloc, Ek,  nn[1], Δασ[1], Aασ_above[1], Aασ_below[1], αασ[1][1], βασ[1][1],  G12ασ[1], Zασ[1]]
n2file="./data_one_band_half/result_n2.dat"
# data=[U,  Etotal,  Eloc, Ek,  nn[1],Zασ[1]]

dmftfile="./dmft_data/one_band_inf_u_d_half.dat"

# the entries for n3 file
set xrange [0:10]
set yrange [0:0.25]
set xlabel "U/t"
set ylabel "d"
plot dmftfile  u 1:2 t "DMFT(NRG)", n2file u 1:5 w l t "N=2",n3file u 1:5 w l t "N=3",

set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_half_d_U.png"
replot
unset terminal


dmftfilez="./dmft_data/one_band_inf_u_Z_half.dat"
set xrange [0:10]
set yrange [0:1]
set xlabel "U/t"
set ylabel "Z"
plot dmftfilez  u 1:2 t "DMFT(NRG)", n2file u 1:6 w l t "N=2",n3file u 1:12 w l t "N=3",

set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_half_Z_U.png"
replot
unset terminal
