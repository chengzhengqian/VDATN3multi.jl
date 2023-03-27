function cal_g011Incr(n,delta,r,c)
(c*(-4.0*n*Power(r,2.0)*(-1.0 + Power(r,2.0)) + 4.0*Power(n,2.0)*Power(r,2.0)*(-1.0 + c + (1.0 + c)*Power(r,2.0)) + delta*(-1.0 + Power(r,2.0))*(4.0*Power(r,2.0) + delta*(-1.0 + Power(r,2.0))*(-1.0 + c + (1.0 + c)*Power(r,2.0)))))/(8.0*c*delta*Power(r,2.0)*(-1.0 + Power(r,2.0)) + Power(delta,2.0)*Power(-1.0 + Power(r,2.0),2.0)*(Power(-1.0 + c,2.0) + Power(1.0 + c,2.0)*Power(r,2.0)) + 4.0*(Power(r,2.0) + Power(r,4.0)) + 4.0*n*Power(r,2.0)*(2.0*(-1.0 + c - (1.0 + c)*Power(r,2.0)) + n*(Power(-1.0 + c,2.0) + Power(1.0 + c,2.0)*Power(r,2.0))))
end
