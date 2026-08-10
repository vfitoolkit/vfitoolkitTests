function F=ReturnFn_nod_z_noe_nosemiz_with2A(aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi)
% Simplified Kitao (2008), no taxes. Two endogenous states:
%   a : assets,  e : occupation (0=worker, 1=entrepreneur)
% eprime is the chosen next-period occupation (free to switch).
% z: eta (labor productivity), theta (entrepreneurial ability).
% (This is BIHAinequality2_ReturnFn with tau_I=tau_c=0, and the leverage
%  parameter renamed from 'd' to 'leverage'.)

F=-Inf;

if e==0 % Worker
    c=(1+r)*a+w*eta-aprime;
    if c>0 && aprime>=0
        F=(c^(1-sigma))/(1-sigma);
    end

elseif e==1 % Entrepreneur (static production problem, solved analytically)
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
    n=((upsilon2*theta*k^upsilon1)/w)^(1/(1-upsilon2));
    output=theta*(k^upsilon1)*(n^upsilon2);
    rbar=r;
    if k>a
        rbar=r+phi; % cost of borrowing
    end
    I=output-delta*k-rbar*(k-a)-w*(n-eta)*(n>eta); % own labor eta is 'free', only pay hired (n-eta)
    c=I+a-aprime; % no taxes, so profit=I+a and c=profit-aprime
    if c>0 && aprime>=0
        F=(c^(1-sigma))/(1-sigma);
    end
end

end
