set terminal qt size 800,600  font ",20"  linewidth 1

n3file(t1)="./two_band_half_t1_t2/two_band_t1_".t1."_t2_1.0_half_result.dat"
#     data=[U,  Etotal,  Eloc, Ek,  nn..., Δασ[[1,3]]..., Aασ_above[[1,3]]..., Aασ_below[[1,3]]..., G12ασ[[1,3]]..., Zασ[[1,3]]...]

# unset label
set xrange [0:7]
set yrange [0:1]
set xlabel "U/t_2"
set ylabel "Z_{/Symbol a}"
set key right top

unset label
set label "t1/t2=0.1" font ",8" at 0.1,0.1
set label "t1/t2=0.2" font ",8" at 2.3,0.4
set label "t1/t2=0.1" font ",8" at 1.4,0.93
set label "t1/t2=0.2" font ",8" at 1.6,0.66



t1s="0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2"
plot NaN lc "blue" t "Z_1",NaN lc "red" t "Z_2", for [t1 in t1s] n3file(t1) u 1:19 w l lc  "blue"   t "", for [t1 in t1s] n3file(t1) u 1:20 w l  lc "red" t  ""

set terminal pngcairo size 800,600 linewidth 1 font ",20"
set output "./figures/two_band_half_inf_Z_U_t1_t2_0.1_0.2.png"
replot
unset terminal


unset label
set xrange [0:7]
set yrange [0:1]
set xlabel "U/t_2"
set ylabel "Z_{/Symbol a}"
set key at 7, 0.6
set label "t1/t2=0.2" font ",20" at 1,0.2
set label "t1/t2=0.3" font ",20" at 2.0,0.8
t1s="0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3"
plot  NaN lc "blue" t "Z_1",NaN lc "red" t "Z_2",  for [t1 in t1s] n3file(t1) u 1:19 w l lc  "blue"   t "", for [t1 in t1s] n3file(t1) u 1:20 w l  lc "red" t  ""

set terminal pngcairo size 800,600 linewidth 1 font ",20"
set output "./figures/two_band_half_inf_Z_U_t1_t2_0.2_0.3.png"
replot
unset terminal



set xrange [0:10]
set yrange [0:1]
set xlabel "U/t_2"
set ylabel "Z_{/Symbol a}"
set key right top

t1s="0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
unset label
set label "t1/t2=0.2" font ",20" at 0.5,0.2
set label "t1/t2=1.0" font ",20" at 4.1,0.7
plot NaN lc "blue" t "Z_1",NaN lc "red" t "Z_2", for [t1 in t1s] n3file(t1) u 1:19 w l lc  "blue"   t "", for [t1 in t1s] n3file(t1) u 1:20 w l  lc "red" t  ""


set terminal pngcairo size 800,600 linewidth 1 font ",20"
set output "./figures/two_band_half_inf_Z_U_t1_t2_0.2_1.0.png"
replot
unset terminal


# examine the double occupancy

n3file(t1)="./two_band_half_t1_t2/two_band_t1_".t1."_t2_1.0_half_result.dat"
#     data=[U,  Etotal,  Eloc, Ek,  nn..., Δασ[[1,3]]..., Aασ_above[[1,3]]..., Aασ_below[[1,3]]..., G12ασ[[1,3]]..., Zασ[[1,3]]...]
# (1, 2, 1.0) (3, 4, 1.0) (1, 4, 1.0) (2, 3, 1.0) (1, 3, 1.0) (2, 4, 1.0), order of interaction pairs

unset label
set xrange [0:10]
set yrange [0:1]
set xlabel "U/t_2"
set ylabel "<n_i*n_j>"

t1s="0.9"
set xrange [0:10]
set yrange [0:0.25]
set key left bottom
set label "t1/t2=0.9" font ",20" at 1.5,0.15
plot NaN lc "red" t "n_{1up}n_{1dn}", NaN lc "green" t "n_{2up}n_{2dn}",NaN lc "blue" t "n_{1{/Symbol s}}n_{2{/Symbol s}}",  for [t1 in t1s] n3file(t1) u 1:5 w l lc "red" t "",  for [t1 in t1s] n3file(t1) u 1:6 w l lc "green" t "" , for [t1 in t1s] n3file(t1) u 1:7 w l lc "blue" t ""


set terminal pngcairo size 800,600 linewidth 1 font ",20"
set output "./figures/two_band_half_inf_nn_U_t1_t2_".t1s."_1.0.png"
replot
unset terminal


t1s="0.5"
t1s="0.3"
set xrange [0:10]
set yrange [0:0.25]
set key left bottom
unset label
set label "t1/t2=".t1s  font ",20" at 6,0.15
plot NaN lc "red" t "n_{1up}n_{1dn}", NaN lc "green" t "n_{2up}n_{2dn}",NaN lc "blue" t "n_{1{/Symbol s}}n_{2{/Symbol s}}",  for [t1 in t1s] n3file(t1) u 1:5 w l lc "red" t "",  for [t1 in t1s] n3file(t1) u 1:6 w l lc "green" t "" , for [t1 in t1s] n3file(t1) u 1:7 w l lc "blue" t ""


set terminal pngcairo size 800,600 linewidth 1 font ",20"
set output "./figures/two_band_half_inf_nn_U_t1_t2_".t1s."_1.0.png"
replot
unset terminal



t1s="0.25"
t1s="0.2"
t1s="0.15"
t1s="0.1"
set xrange [0:10]
set yrange [0:0.25]
set xlabel "U/t_2"
set ylabel "<n_i*n_j>"
set key at 9.8, 0.1
unset label
set label "t1/t2=".t1s  font ",20" at 6,0.15
plot NaN lc "red" t "n_{1up}n_{1dn}", NaN lc "green" t "n_{2up}n_{2dn}",NaN lc "blue" t "n_{1{/Symbol s}}n_{2{/Symbol s}}",  for [t1 in t1s] n3file(t1) u 1:5 w l lc "red" t "",  for [t1 in t1s] n3file(t1) u 1:6 w l lc "green" t "" , for [t1 in t1s] n3file(t1) u 1:7 w l lc "blue" t ""


set terminal pngcairo size 800,600 linewidth 1 font ",20"
set output "./figures/two_band_half_inf_nn_U_t1_t2_".t1s."_1.0.png"
replot
unset terminal

