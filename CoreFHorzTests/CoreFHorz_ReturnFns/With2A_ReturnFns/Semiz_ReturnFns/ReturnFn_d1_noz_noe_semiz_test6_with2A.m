function F=ReturnFn_d1_noz_noe_semiz_test6_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2)
% Cross-test 6 helper. d1 is in the signature but does NOT enter the
% return value, so the with-d1 model should give the same V (and any d1
% choice is optimal) as the nod1 model with the same other inputs.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*semiz+uempbenefit*(1-semiz)-a1prime+a2-a2prime;
else
    c=(1+r)*a1+pension-a1prime+a2-a2prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d2;
    if agej<Jr
        F=F+phi1*a2^phi2;
    end
end

end
