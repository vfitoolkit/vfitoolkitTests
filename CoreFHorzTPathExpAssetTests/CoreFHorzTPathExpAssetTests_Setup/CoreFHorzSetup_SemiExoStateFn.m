function prob=CoreFHorzSetup_SemiExoStateFn(n,nprime,dsemiz,probfindjob,problosejob)

prob=-1; % placeholder

% dsemiz=1 makes states very highly persistent
if dsemiz==1
    probfindjob=0.1*probfindjob;
    problosejob=0.1*problosejob;
end

if n==0
    if nprime==1
        prob=probfindjob;
    elseif nprime==0
        prob=1-probfindjob;
    end
elseif n==1
    if nprime==1
        prob=1-problosejob;
    elseif nprime==0
        prob=problosejob;
    end
end

end