function F=ReturnFn_nod1_z_e_noa1_semiz_with2A2(d2,d3,a2_1,a2_2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate,uempbenefit,searcheffortcost)

F=-Inf;

if agej<Jr
    c=w*kappa_j*d2*a2_1*z*e*semiz + uempbenefit*(1-semiz);
else
    c=pension+pensionrate*a2_2;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end
