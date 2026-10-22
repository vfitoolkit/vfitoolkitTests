function F=ReturnFn_nod_z_noe_InfNaNcap(aprime,a,z,r,w,sigma,acap,zcap)
% ReturnFn_nod_z_noe_nosemiz plus a cap on the STATE: for z above zcap, any a above acap is -Inf whatever
% aprime is, so V(a,z)=-Inf there while V(a,z) at the other z is finite. This puts a dead region ABOVE
% finite states (the natural borrowing limit only puts one below), which is what makes the grid
% interpolation layer's Howards step put a zero weight on a -Inf UPPER node. See doPart(11) of
% CoreInfHorzVFIAlgoTests and InfHorzGI_InfNaN_proposal.md at the toolkit root.

F=-Inf;

if a<=acap || z<=zcap
    c=(1+r)*a+w*z-aprime;
    if c>0
        F=(c^(1-sigma)-1)/(1-sigma);
    end
end

end
