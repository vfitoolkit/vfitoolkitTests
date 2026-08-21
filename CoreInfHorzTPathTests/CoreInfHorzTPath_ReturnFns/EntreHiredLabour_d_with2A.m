function nhired=EntreHiredLabour_d_with2A(l,aprime,eprime,a,e,eta,theta,r,w,delta,upsilon1,upsilon2,leverage,phi)
% Note: l, aprime and eprime are unused here, but the toolkit convention is that the function
% called from a FnsToEvaluate anonymous fn carries the state variables through anyway.
% Labour the entrepreneur HIRES IN from the market, for the 2A (Kitao) model.
% Zero for workers (e==0).
%
% The entrepreneur's total labour input is n. Their own eta is used first and is free (the
% ReturnFn pays only w*(n-eta) when n>eta), so the amount taken out of the market is max(n-eta,0).
% This is why the corporate sector's labour is (workers' labour supply) minus (this), rather than
% total labour supply minus n.
%
% As with EntreCapital_with2A, the k and n blocks REPLICATE those inside ReturnFn_*_with2A and
% must be kept in step with them.

nhired=0;
if e==1
    onegg=1-upsilon1-upsilon2;
    k_unc=(upsilon1/(r+phi+delta))^((1-upsilon1)/onegg) * (upsilon2/w)^(upsilon2/onegg) * theta^(1/onegg);
    if k_unc>(1+leverage)*a
        k=(1+leverage)*a;
    else
        k=k_unc;
    end
    if k<a
        k_r=(upsilon1/(r+delta))^((1-upsilon1)/onegg) * (upsilon2/w)^(upsilon2/onegg) * theta^(1/onegg);
        if k_r>a
            k=a;
        else
            k=k_r;
        end
    end
    n=((upsilon2*theta*k^upsilon1)/w)^(1/(1-upsilon2));
    if n>eta
        nhired=n-eta;
    end
end

end
