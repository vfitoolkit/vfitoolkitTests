function F=ReturnFn_d_z_e_InfNaNcap(d,aprime,a,z,e,r,w,sigma,eta,varphi,acap,zcap)
% ReturnFn_d_z_e_nosemiz plus the cap on the state (a above acap at z above zcap is -Inf), see
% ReturnFn_nod_z_noe_InfNaNcap. The cap does not depend on e.

F=-Inf;

if a<=acap || z<=zcap
    c=(1+r)*a+w*d*z*e-aprime;
    if c>0 && d<1
        F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d)^(1-eta)-1)/(1-eta);
    end
end

end
