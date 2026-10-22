function F=ReturnFn_d1_z_noe_InfNaNcap_semiz(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost,acap,zcap)
% ReturnFn_d1_z_noe_semiz plus the cap on the state (a above acap at z above zcap is -Inf), see ReturnFn_nod_z_noe_InfNaNcap.

F=-Inf;

if a<=acap || z<=zcap
    c=(1+r)*a+w*d1*z*semiz+uempbenefit*(1-semiz)-aprime;
    if c>0 && d1<1
        F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta)-searcheffortcost*d2;
    end
end

end
