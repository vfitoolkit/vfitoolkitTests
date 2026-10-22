function F=ReturnFn_nod_z_e_InfNaNcap(aprime,a,z,e,r,w,sigma,acap,zcap)
% ReturnFn_nod_z_e_nosemiz plus the cap on the state (a above acap at z above zcap is -Inf), see
% ReturnFn_nod_z_noe_InfNaNcap. The cap does not depend on e.

F=-Inf;

if a<=acap || z<=zcap
    c=(1+r)*a+w*z*e-aprime;
    if c>0
        F=(c^(1-sigma)-1)/(1-sigma);
    end
end

end
