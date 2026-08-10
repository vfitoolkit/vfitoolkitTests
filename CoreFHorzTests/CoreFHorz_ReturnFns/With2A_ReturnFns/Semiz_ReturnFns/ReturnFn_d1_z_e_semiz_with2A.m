function F=ReturnFn_d1_z_e_semiz_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*semiz*z*e+uempbenefit*(1-semiz)-a1prime+a2-a2prime;
else
    c=(1+r)*a1+pension-a1prime+a2-a2prime;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta)-searcheffortcost*d2;
    if agej<Jr
        F=F+phi1*a2^phi2; % give utility from a2 during working ages
    end
end


end
