function F=ReturnFn_d_z_noe_InfNaNcap_with2A(d,a1prime,a2prime,a1,a2,z,r,r2,w,sigma,eta,varphi,acap,zcap)
% ReturnFn_d_z_noe_nosemiz_with2A plus the cap on the first endogenous state (a1 above acap at z
% above zcap is -Inf), see ReturnFn_nod_z_noe_InfNaNcap_with2A.

F=-Inf;

if a1<=acap || z<=zcap
    c=(1+r)*a1+(1+r2)*a2+w*d*z-a1prime-a2prime;
    if c>0 && d<1
        F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d)^(1-eta)-1)/(1-eta);
    end
end

end
