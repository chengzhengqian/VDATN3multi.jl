function cal_dnAdab(a,b,e)
[(Power(b,2.0)/(2.0*Power(Power(b,2.0) + Power(a - e,2.0),1.5))),  ((b*(-a + e))/(2.0*Power(Power(b,2.0) + Power(a - e,2.0),1.5))),  ((b*(-a + e))/(2.0*Power(Power(b,2.0) + Power(a - e,2.0),1.5))),  (Power(a - e,2.0)/(2.0*Power(Power(b,2.0) + Power(a - e,2.0),1.5))),  ]
end
