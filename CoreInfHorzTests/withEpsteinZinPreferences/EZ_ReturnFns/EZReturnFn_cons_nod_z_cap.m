function F=EZReturnFn_cons_nod_z_cap(aprime,a,z,r,w,acap,zcap)
% EZReturnFn_cons_nod_z plus a cap on the STATE: for z above zcap, any a above acap is -Inf whatever
% aprime is, so V(a,z)=-Inf there while V(a,z) at the other z is finite. This puts a dead region ABOVE
% finite states (the natural borrowing limit only puts one below), which is what makes the grid
% interpolation layer's Howards step put a zero weight on a -Inf UPPER node. Mirrors
% ReturnFn_nod_z_noe_InfNaNcap of CoreInfHorzVFIAlgoTests (InfHorzGI_InfNaN_proposal.md).

F=-Inf;

if a<=acap || z<=zcap
    c=(1+r)*a+w*z-aprime;
    if c>0
        F=c;
    end
end

end
