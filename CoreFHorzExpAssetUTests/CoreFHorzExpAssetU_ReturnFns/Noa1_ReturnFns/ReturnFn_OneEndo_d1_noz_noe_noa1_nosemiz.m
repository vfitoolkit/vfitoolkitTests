function F=ReturnFn_OneEndo_d1_noz_noe_noa1_nosemiz(d1,aprime,a,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension)
% Standard 1-endo Case1 (with d1). Used as side A of CrossTests3 noa1 d1 version.

F=-Inf;

if agej<Jr
    c=w*kappa_j*d1*aprime*a;
else
    c=pension;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta);
end


end
