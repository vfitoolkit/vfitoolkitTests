function prob=CoreFHorzTwoEndoSetup_SemiExoStateFn_TwoMarkovs(n,nprime,dsemiz, ...
    p_odd_00, p_odd_01, p_odd_10, p_odd_11, p_even_00, p_even_01, p_even_10, p_even_11)
% Cross-test helper. Two 2x2 markov matrices selected by dsemiz.
% semiz_grid = [0;1]. dsemiz=1 -> pi_odd, dsemiz=2 -> pi_even.
% p_X_YZ is the (Y,Z) entry (Y=current state, Z=next state, with 0/1 = grid value).

prob=0;
if dsemiz==1
    if n==0 && nprime==0, prob=p_odd_00; end
    if n==0 && nprime==1, prob=p_odd_01; end
    if n==1 && nprime==0, prob=p_odd_10; end
    if n==1 && nprime==1, prob=p_odd_11; end
elseif dsemiz==2
    if n==0 && nprime==0, prob=p_even_00; end
    if n==0 && nprime==1, prob=p_even_01; end
    if n==1 && nprime==0, prob=p_even_10; end
    if n==1 && nprime==1, prob=p_even_11; end
end

end
