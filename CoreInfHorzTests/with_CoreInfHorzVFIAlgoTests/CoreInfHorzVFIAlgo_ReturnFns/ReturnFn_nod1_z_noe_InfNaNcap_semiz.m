function F=ReturnFn_nod1_z_noe_InfNaNcap_semiz(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost,acap,zcap)
% ReturnFn_nod1_z_noe_semiz plus the cap on the state (a above acap at z above zcap is -Inf, whatever d2, aprime or
% semiz), see ReturnFn_nod_z_noe_InfNaNcap.

F=-Inf;

if a<=acap || z<=zcap
    c=(1+r)*a+w*z*semiz+uempbenefit*(1-semiz)-aprime;
    if c>0
        F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d2;
    end
end

end
