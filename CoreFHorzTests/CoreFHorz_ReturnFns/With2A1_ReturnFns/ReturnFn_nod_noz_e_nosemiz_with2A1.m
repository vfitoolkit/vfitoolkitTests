function F=ReturnFn_nod_noz_e_nosemiz(a1prime,a2prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*e-a1prime+a2-a2prime;
else
    c=(1+r)*a1+pension-a1prime+a2-a2prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
    if agej<Jr
        F=F+phi1*a2^phi2; % give utility from a2, but only during working ages (saves me doing a decent job of setting up model, as households want some amount of a2, but sell it all when they retire)
    end
end


end