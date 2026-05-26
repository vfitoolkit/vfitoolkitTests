function F=ReturnFn_OneEndo_nod1_noz_noe_noa1_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension)
% Standard 1-endo Case1 (no d). Used as side A of CrossTests3 noa1: aprime plays the role of the
% experience asset's "next-period chosen value" (the ExpAsset side has aprimeFn=@(d2,a2) d2, so d2=aprime).

F=-Inf;

if agej<Jr
    c=w*kappa_j*aprime*a;
else
    c=pension;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
