function F=ReturnFn_nod1_noz_noe_semiz_test2_with2A(d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2)
% Cross-test 2 helper. Same consumption as the z-only model. Forces d2=1
% at odd ages and d2=2 at even ages by returning -Inf otherwise. With this
% forcing, pi_semiz(:,:,d2,jj) selects pi_odd at odd ages and pi_even at
% even ages, matching the alternating pi_z_J in the markov-z comparison.

F=-Inf;

% Force d2 by age parity before computing F
if mod(agej,2)==1 && d2~=1
    return
end
if mod(agej,2)==0 && d2~=2
    return
end

if agej<Jr
    c=(1+r)*a1+w*kappa_j*semiz-a1prime+a2-a2prime;
else
    c=(1+r)*a1+pension-a1prime+a2-a2prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
    if agej<Jr
        F=F+phi1*a2^phi2;
    end
end

end
