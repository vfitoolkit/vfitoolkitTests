function prob=CoreFHorzSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2)

prob=-1; % placeholder

% Notice that dsemiz does nothing, hence this is actually a markov,
% disguised as a semi-exogenous state.
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