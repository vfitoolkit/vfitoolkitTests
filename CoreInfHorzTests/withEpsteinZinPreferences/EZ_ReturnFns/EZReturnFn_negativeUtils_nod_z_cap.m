function F=EZReturnFn_negativeUtils_nod_z_cap(aprime,a,z,r,w,ezsigma,acap,zcap)
% EZReturnFn_negativeUtils_nod_z plus a cap on the STATE (see EZReturnFn_cons_nod_z_cap).

F=-Inf;

if a<=acap || z<=zcap
    c=(1+r)*a+w*z-aprime;
    if c>0
        F=(c^(1-ezsigma))/(1-ezsigma);
    end
end

end
