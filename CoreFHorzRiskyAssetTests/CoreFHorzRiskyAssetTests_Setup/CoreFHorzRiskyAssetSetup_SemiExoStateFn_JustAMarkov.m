function prob=CoreFHorzRiskyAssetSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2)
% Semi-exogenous state whose transitions do NOT depend on the decision dsemiz,
% so it is really just a markov disguised as a semi-exogenous state. Used by the
% semiz-as-z cross test: with this SemiExoStateFn (and searcheffortcost=0 so dsemiz
% is inert) the semiz model reproduces the z-markov model exactly.
% States take grid values z1,z2. Transition matrix:
%   from z1: [1-probfindjob (stay z1), probfindjob (to z2)]
%   from z2: [problosejob (to z1),     1-problosejob (stay z2)]

prob=-1; % placeholder

if n==z1
    if nprime==z2
        prob=probfindjob;
    elseif nprime==z1
        prob=1-probfindjob;
    end
elseif n==z2
    if nprime==z2
        prob=1-problosejob;
    elseif nprime==z1
        prob=problosejob;
    end
end

end
