function F=EZReturnFn_positiveUtils_d_z(d,aprime,a,z,r,w,varphi,ezsigma)
% Utility-units, POSITIVE-valued utility fn for Epstein-Zin
% (vfoptions.EZutils=1, vfoptions.EZpositiveutility=1).
% The (1+x) shift keeps F strictly positive for all x>0 [a plain
% (x^(1-ezsigma)-1)/(1-ezsigma) would change sign at x=1].
% InfHorz form of EZReturnFn_positiveUtils_d_z_noe_nosemiz (no ages, no retirement).

F=-Inf;

c=(1+r)*a+w*d*z-aprime;

if c>0 && d<1
    x=(c^varphi)*((1-d)^(1-varphi));
    F=((1+x)^(1-ezsigma)-1)/(1-ezsigma);
end


end
