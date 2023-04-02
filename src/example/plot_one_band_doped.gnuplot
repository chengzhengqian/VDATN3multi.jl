set terminal qt size 800,600  font ",20"

n3file(n)="./data_one_band_doped/result_n3_ntotal_".n.".dat"
n2file(n)="./data_one_band_doped/result_n2_ntotal_".n.".dat"
dmftfile(n)="./dmft_data/one_band_inf_u_d_n_".n.".dat"




set xrange [0:10]
set yrange [0:0.19]
set xlabel "U/t"
set ylabel "d"
set key at graph 0.9, graph 0.95
ns="0.4 0.6 0.8 0.9"
plot for [n in ns ] n2file(n) u 1:5 w l t "" lc "green", for [n in ns ] n3file(n) u 1:5 w l t "" lc "red",  for [n in ns ] dmftfile(n)  u 1:2 t ""  pt 5 lc "black", NaN  t "DMFT(NRG)" w p pt 5  lc "black", NaN w l lc "green" t "N=2", NaN w l lc "red" t "N=3"


set label "n=0.4" textcolor "black" at first 0.3, first 0.01
set label "n=0.6" textcolor "black" at first 0.3, first 0.04
set label "n=0.8" textcolor "black" at first 0.3, first 0.09
set label "n=0.9" textcolor "black" at first 2.4, first 0.14
# unset label

set terminal pngcairo size 800,600 linewidth 3 font ",20"
set output "./figures/one_band_doped_d_U.png"
replot
unset terminal

