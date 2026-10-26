function F=EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma)
% Utility-units, POSITIVE-valued utility fn for Epstein-Zin
% (vfoptions.EZutils=1, vfoptions.EZpositiveutility=1).
% The (1+x) shift keeps F strictly positive for all x>0 [a plain
% (x^(1-ezsigma)-1)/(1-ezsigma) would change sign at x=1].
% InfHorz form of EZReturnFn_positiveUtils_nod_z_noe_nosemiz (no ages, no retirement).

F=-Inf;

c=(1+r)*a+w*z-aprime;

if c>0
    x=c;
    F=((1+x)^(1-ezsigma)-1)/(1-ezsigma);
end


end
