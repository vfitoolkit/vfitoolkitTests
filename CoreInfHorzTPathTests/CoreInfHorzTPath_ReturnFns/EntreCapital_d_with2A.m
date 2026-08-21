function k=EntreCapital_d_with2A(l,aprime,eprime,a,e,theta,r,w,delta,upsilon1,upsilon2,leverage,phi)
% Note: l, aprime and eprime are unused here, but the toolkit convention is that the function
% called from a FnsToEvaluate anonymous fn carries the state variables through anyway.
% Capital used in the entrepreneur's own business, for the 2A (Kitao) model.
% Zero for workers (e==0).
%
% This REPLICATES the k block inside ReturnFn_*_with2A. The duplication is unavoidable: the
% household problem needs k to form the return, and the general eqm needs k aggregated to split
% total assets into non-corporate and corporate capital. The Kitao/Bruggemann replications do the
% same. If the k block in the ReturnFn is ever edited, this must be edited to match.

k=0;
if e==1
    onegg=1-upsilon1-upsilon2;
    % Unconstrained capital if borrowing (rate r+phi)
    k_unc=(upsilon1/(r+phi+delta))^((1-upsilon1)/onegg) * (upsilon2/w)^(upsilon2/onegg) * theta^(1/onegg);
    if k_unc>(1+leverage)*a % collateral constraint
        k=(1+leverage)*a;
    else
        k=k_unc;
    end
    if k<a % if not borrowing, use the saving rate r
        k_r=(upsilon1/(r+delta))^((1-upsilon1)/onegg) * (upsilon2/w)^(upsilon2/onegg) * theta^(1/onegg);
        if k_r>a
            k=a;
        else
            k=k_r;
        end
    end
end

end
