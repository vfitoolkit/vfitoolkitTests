function aprime=aprimeFn_CoreTestRiskyAsset(riskyshare,savings,u,r)
% Portfolio-choice riskyasset.
% Inputs are (d2,d3,u,...) under vfoptions.refine_d — d2=riskyshare, d3=savings.
% Decomposes savings into safe (return r) and risky (return r+u) parts.

if savings>0
    aprime=(1+r)*(1-riskyshare)*savings + (1+r+u)*riskyshare*savings;
else
    % negative savings ≡ mortgage; cannot be invested risky, force safe
    aprime=(1+r)*savings;
end

end
