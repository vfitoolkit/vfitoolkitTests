function jequaloneDist=EstimFHorz_jequaloneDistParamFn(a_grid,z_grid,n_a,n_z,initdecay)
% jequaloneDist-as-a-function WITH a parameter (EstimFHorzTests part 9): geometric mass with decay initdecay over the first SIXTY
% asset gridpoints, all at the middle z point. In part 9 initdecay is DERIVED by the ParametrizeParamsFn from the estimated
% loginitdecay, pinning down that the objective function applies ParametrizeParamsFn before evaluating jequaloneDist.
% Sixty points because the bank's a_grid is 5*linspace(0,1,101).^3: its first ten points all sit within 0.004 of zero, so a decay
% parameter spread over only those is economically invisible (the first version of this test was unidentified for exactly that
% reason -- the optimiser faced a flat objective). Over sixty points (up to ~1.0) the decay moves real wealth.
% (a_grid and z_grid are unused, but the toolkit's jequaloneDist-function convention passes them first.)

w=initdecay.^(0:59)';
w=w/sum(w);
jequaloneDist=zeros(n_a,prod(n_z));
jequaloneDist(1:60,ceil(prod(n_z)/2))=w;

end
