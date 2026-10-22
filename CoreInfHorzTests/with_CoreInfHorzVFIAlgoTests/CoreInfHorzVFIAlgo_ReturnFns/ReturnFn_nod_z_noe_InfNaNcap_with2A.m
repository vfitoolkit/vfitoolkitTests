function F=ReturnFn_nod_z_noe_InfNaNcap_with2A(a1prime,a2prime,a1,a2,z,r,r2,w,sigma,acap,zcap)
% ReturnFn_nod_z_noe_nosemiz_with2A plus the cap on the FIRST endogenous state (a1 above acap at z
% above zcap is -Inf), see ReturnFn_nod_z_noe_InfNaNcap. The grid interpolation layer is on a1, so
% that is where the dead region has to sit; a2 is free.

F=-Inf;

if a1<=acap || z<=zcap
    c=(1+r)*a1+(1+r2)*a2+w*z-a1prime-a2prime;
    if c>0
        F=(c^(1-sigma)-1)/(1-sigma);
    end
end

end
