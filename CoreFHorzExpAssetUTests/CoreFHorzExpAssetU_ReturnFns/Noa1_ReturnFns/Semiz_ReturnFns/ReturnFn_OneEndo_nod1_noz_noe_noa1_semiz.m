function F=ReturnFn_OneEndo_nod1_noz_noe_noa1_semiz(d3,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost,transfer)
% NOTE (experienceassetu noa1 only): working-age consumption is earnings alone, and a is the
% experience asset, which is not consumable. At the first grid point a=0 earnings are exactly
% 0, so without 'transfer' F=-Inf there; the experience-asset interpolation then puts positive
% weight on that point (weight 1 when a2'<=a2_grid(1)), so the -Inf propagated to every state
% and backwards through every age, making whole blocks of this tier's checks vacuous.
% A transfer, not a max(c,floor): a floor makes F constant below it and manufactures exact
% ties in the argmax, which this project has had enough trouble with already.
% Standard 1-endo Case1 with semiz (no d2 since aprime is the choice). Used as side A of
% CrossTests3 noa1 semiz: aprime plays the role of d2 in the ExpAsset side (where aprimeFn=@(d2,a2) d2).

F=-Inf;

if agej<Jr
    c=w*kappa_j*aprime*a*semiz+uempbenefit*(1-semiz)+transfer; % transfer keeps c>0 even at a=0, where working-age earnings are exactly 0 because a is the (unconsumable) experience asset: see the header note
else
    c=pension;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end
