function F=ReturnFn_ExpAssetze_d1_z_e_semiz_noa1(d1,d2,d3,a2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost)
% noa1 version of ReturnFn_ExpAssetze_d1_z_e_semiz: no a1prime/a1 args (the
% experience asset a2 is the only endogenous state, so no savings decision)

F=-Inf;

if agej<Jr
    c=w*kappa_j*d1*d2*a2*semiz*z*e+uempbenefit*(1-semiz);
else
    c=pension;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta)-searcheffortcost*d3;
end


end
