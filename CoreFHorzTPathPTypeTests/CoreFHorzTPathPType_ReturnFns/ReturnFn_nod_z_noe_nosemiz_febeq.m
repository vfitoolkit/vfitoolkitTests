function F=ReturnFn_nod_z_noe_nosemiz_febeq(aprime,a,z,r,w,fe,kappa_j,sigma,agej,Jr,pension,beq)
% As ReturnFn_nod_z_noe_nosemiz, but earnings carry a fixed effect fe and
% everyone receives a lump-sum bequest transfer beq.

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*fe*z+beq-aprime;
else
    c=(1+r)*a+pension+beq-aprime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
