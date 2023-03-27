function cal_glocInnA(n,delta,Abl,Aab,sloc11,sloc12)
[(n)  ((delta*sloc11)/sloc12)  ((Aab*sloc11)/Sqrt(sloc12*(Power(sloc11,2.0) + Power(sloc12,2.0))))  ;
(-((delta*sloc11)/sloc12))  (n)  (Abl*Sqrt(1.0/sloc12) + Aab/Sqrt(Power(sloc11,2.0)/sloc12 + sloc12))  ;
(-(Abl*sloc11*Sqrt(1.0/sloc12)))  (-(Abl*Sqrt(sloc12)) - Aab*Sqrt(Power(sloc11,2.0)/sloc12 + sloc12))  (n)  ;
]
end
