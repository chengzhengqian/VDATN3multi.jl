set terminal qt size 800,600  font ",20"  linewidth 1

n3file(t1,J_U)="./two_band_half_t1_t2_with_J/two_band_t1_".t1."_t2_1.0_J_U_".J_U."_half_result.dat"



#     data=[U,  Etotal,  Eloc, Ek,  nn..., Δασ[[1,3]]..., Aασ_above[[1,3]]..., Aασ_below[[1,3]]..., G12ασ[[1,3]]..., Zασ[[1,3]]...]

set xrange [0:6]
set yrange [0:1]
set xlabel "U/t_2"
set ylabel "Z_{/Symbol a}"
set key right top

t1s="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
# J_U="0.05"
# J_U="0.1"
J_U="0.25"
set key at 5.5 ,0.95
unset label
# J_U=0.05
# set label "t1/t2=1.0" font ",10" at 4.0,0.65
# set label "t1/t2=0.1" font ",10" at 1.0,0.1
# J_U=0.1
set label "t1/t2=1.0" font ",10" at 3.5,0.65
set label "t1/t2=0.1" font ",10" at 0.8,0.1
# J_U=0.25
set xrange [0:4]
set yrange [0:1]
unset label
set key at 3.9 ,0.95
set label "t1/t2=1.0" font ",10" at 3.0,0.65
set label "t1/t2=0.1" font ",10" at 0.6,0.1
set label  "J/U=".J_U font ",20" at 2.0,0.90


plot NaN lc "blue" t "Z_1",NaN lc "red" t "Z_2", for [t1 in t1s] n3file(t1,J_U) u 1:19 w l lc  "blue"   t "", for [t1 in t1s] n3file(t1,J_U) u 1:20 w l  lc "red" t  ""

set terminal pngcairo size 800,600 linewidth 1 font ",20"
set output "./figures/two_band_half_inf_Z_U_t1_t2_0.1_1.0_J_U_".J_U.".png"
replot
unset terminal

unset label
set xrange [0:4]
set yrange [0:0.5]
set xlabel "U/t_2"
set ylabel "<n_i*n_j>"

unset label
t1s="0.9"
t1s="0.5"
t1s="0.1"
J_U="0.25"
set xrange [0:4]
set yrange [0:0.5]
set key at 0.5 ,0.3
# t1s="0.9"
set label "J/U=".J_U  font ",20" at 0.5,0.15
set label "t1/t2=".t1s font ",20" at 0.5,0.1
# t1s="0.5"
set label "J/U=".J_U  font ",20" at 0.5,0.1
set label "t1/t2=".t1s font ",20" at 0.5,0.05
# t1s="0.1"
set key at 0.4 ,0.33
set label "J/U=".J_U  font ",20" at 0.5,0.15
set label "t1/t2=".t1s font ",20" at 0.5,0.1
# (1, 2, 1.0) (3, 4, 1.0) (1, 4, 1.0) (2, 3, 1.0) (1, 3, 1.0) (2, 4, 1.0), order of interaction pairs

plot NaN lc "red" t "n_{1up}n_{1dn}", NaN lc "green" t "n_{2up}n_{2dn}",NaN lc "blue" t "n_{1up}n_{2dn}", NaN lc "black" t "n_{1up}n_{2up}",  for [t1 in t1s] n3file(t1,J_U) u 1:5 w l lc "red" t "",  for [t1 in t1s] n3file(t1,J_U) u 1:6 w l lc "green" t "" , for [t1 in t1s] n3file(t1,J_U) u 1:7 w l lc "blue" t "", for [t1 in t1s] n3file(t1,J_U) u 1:9 w l lc "black" t ""

set terminal pngcairo size 800,600 linewidth 1 font ",20"
set output "./figures/two_band_half_inf_nn_U_t1_t2_".t1s."_1.0_J_U_".J_U.".png"
replot
unset terminal
