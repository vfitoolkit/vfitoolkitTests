function F=ReturnFn_d_z_noe_nosemiz_with2A(l,aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi,xi,sigma2,lbar)
% Simplified Kitao (2008) two-endogenous-state model with Bruggemann (2021)
% endogenous labor supply. No taxes. Two endogenous states:
%   a : assets,  e : occupation (0=worker, 1=entrepreneur)
% Decision variable:
%   l : labor supply (workers choose it; entrepreneurs must supply lbar)
% Separable utility: c^(1-sigma)/(1-sigma) - xi*labsup^(1+sigma2)/(1+sigma2)

F=-Inf;

if e==0 % Worker: chooses labor supply l
    c=(1+r)*a+w*l*eta-aprime;
    if c>0 && aprime>=0
        F=(c^(1-sigma))/(1-sigma) - xi*(l^(1+sigma2))/(1+sigma2);
    end

elseif e==1 % Entrepreneur: must supply l=lbar (all other l are infeasible, giving a unique labor policy)
    if abs(l-lbar)<1e-9
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
        output=theta*(k^upsilon1)*(n^upsilon2);
        rbar=r;
        if k>a
            rbar=r+phi;
        end
        I=output-delta*k-rbar*(k-a)-w*(n-eta)*(n>eta);
        c=I+a-aprime;
        if c>0 && aprime>=0
            F=(c^(1-sigma))/(1-sigma) - xi*(lbar^(1+sigma2))/(1+sigma2);
        end
    end
end

end
