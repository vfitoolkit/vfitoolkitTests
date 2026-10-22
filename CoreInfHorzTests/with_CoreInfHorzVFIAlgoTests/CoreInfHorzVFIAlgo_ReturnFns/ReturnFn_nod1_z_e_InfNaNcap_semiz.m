function F=ReturnFn_nod1_z_e_InfNaNcap_semiz(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost,acap,zcap)
% ReturnFn_nod1_z_e_semiz plus the cap on the state (a above acap at z above zcap is -Inf), see ReturnFn_nod_z_noe_InfNaNcap.

F=-Inf;

if a<=acap || z<=zcap
    c=(1+r)*a+w*z*semiz*e+uempbenefit*(1-semiz)-aprime;
    if c>0
        F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d2;
    end
end

end
