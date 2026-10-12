function F=ReturnFn_d1_noz_e_noa1_nosemiz(d1,d2,a,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,transfer)
% NOTE (experienceassetu noa1 only): working-age consumption is earnings alone, and a is the
% experience asset, which is not consumable. At the first grid point a=0 earnings are exactly
% 0, so without 'transfer' F=-Inf there; the experience-asset interpolation then puts positive
% weight on that point (weight 1 when a2'<=a2_grid(1)), so the -Inf propagated to every state
% and backwards through every age, making whole blocks of this tier's checks vacuous.
% A transfer, not a max(c,floor): a floor makes F constant below it and manufactures exact
% ties in the argmax, which this project has had enough trouble with already.

F=-Inf;

if agej<Jr
    c=w*kappa_j*d1*d2*a*e+transfer; % transfer keeps c>0 even at a=0, where working-age earnings are exactly 0 because a is the (unconsumable) experience asset: see the header note
else
    c=pension;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta);
end


end
