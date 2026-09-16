function DiscP5_runtimetable(outputP5,calib,znums)
% Runtime table for P5 (VAR(1) with normal innovations)
%
% Reported, never asserted. The cost here grows faster than anywhere else in the bank: the state
% count is znum^M and the transition matrix is znum^(2M), so a step in M is not comparable to a
% step in znum. The table is laid out rows-by-znum and columns-by-M to make that visible.

fprintf('\n========== Runtimes for discretizing a VAR(1) with normal innovations ========== \n')
v=version;
try
    g=gpuDevice; gname=g.Name;
catch
    gname='no gpu';
end
fprintf('(median of up to %i timed calls after 1 discarded warm-up; MATLAB %s; GPU: %s) \n',calib.nreps,v,gname)
fprintf('cells with znum^M above the ceiling of %i states are skipped and shown as "skip" \n',calib.Ncap)

nmlist={'discretizeVAR1_Tauchen','discretizeVAR1_FarmerToda'};
fldlist={'VAR1_Tauchen','VAR1_FarmerToda'};
Mlist=outputP5.VAR1_Tauchen.Mlist;
nz=length(znums);
for m_c=1:2
    fprintf('\n--- %s (rows znum, columns M) --- \n',nmlist{m_c})
    fprintf('%8s',' ');
    for c_c=1:length(Mlist)
        fprintf('%12s',['M=',num2str(Mlist(c_c))]);
    end
    fprintf('%14s \n','states at M=2');
    t=outputP5.(fldlist{m_c}).runtime;
    sk=outputP5.(fldlist{m_c}).skipped;
    for c_c=1:nz
        fprintf('%8i',znums(c_c));
        for k_c=1:length(Mlist)
            if sk(c_c,k_c)
                fprintf('%12s','skip');
            else
                fprintf('%12.6f',t(c_c,k_c));
            end
        end
        fprintf('%14i \n',znums(c_c)^2);
    end
end

nskipT=sum(outputP5.VAR1_Tauchen.skipped(:));
fprintf('\n%i of %i cells were skipped for size. \n',nskipT,nz*length(Mlist))
fprintf('Reading the table across a row rather than down a column is the point: going from M=1 to \n')
fprintf('M=2 at fixed znum squares the state count, where going from znum=15 to znum=31 at fixed \n')
fprintf('M=2 multiplies it by about four. The two commands pay for that differently - VAR1_Tauchen \n')
fprintf('calls mvncdf once per row over prod(znum) rectangles, VAR1_FarmerToda solves a maximum \n')
fprintf('entropy problem once per row PER VARIABLE, so its cost carries an extra factor of M. \n')

end
