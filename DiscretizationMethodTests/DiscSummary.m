function DiscSummary(diaryfilename)
% Read the diary back and print one verdict for the whole run.
%
% WHY THIS EXISTS. The bank prints its results in four different formats, and only one of them is
% checkable by eye at a glance:
%    this should be one: 1                       - a flag, fails at 0
%    this should be zero: 2.22e-16               - either a flag or an identity residual
%    this should be below 1e-07: 3.33e-16        - the reader has to do the comparison
%    [MC, pass if |diff|<4*se: 1]                - a different syntax entirely
% With around 2600 checks in a 7600-line diary, "did the run pass" is not a question anyone should
% have to answer by scrolling, and a reader who scans for only one of these forms will report a
% clean run while never having looked at the other three. That happened repeatedly while the bank
% was being built. So the rule lives here, in one place, applied to all four forms.
%
% The diary is closed, read, and reopened in append mode, so this summary lands in the diary itself.
%
% THRESHOLDS. A "should be zero" line is a flag if its value is written as a bare integer and an
% identity residual otherwise; residuals are held to 1e-10, which is the bank's T1 bar. On the run
% of 2026-08-27 the largest residual anywhere was 2.95e-13, so that bar has three orders of
% headroom and will not fire on ordinary floating-point noise.
%
% WHERE A FAILURE IS SAID TO SIT. Every failure is reported under the last ========== banner above
% it, which is always the enclosing subcode. A "--- ... ---" section header is added on top of that
% only when it is genuinely enclosing, and the run of 2026-09-16 is why that qualification exists:
% the wSV_FarmerToda diagonal-convergence check is printed AFTER the calibration loop it summarises,
% so the nearest --- header was "calibration drift (rho=0.9)" from 710 lines earlier, and the
% summary reported the failure as sitting in a calibration it is not part of. Two rules now stop
% that. A --- section is dropped when a ========== banner is crossed, since the banner means a new
% subcode has started, and it is dropped when it is more than seclookback lines back, since a header
% that far up is the last thing that happened rather than the thing this line belongs to.

restol=10^(-10);
seclookback=120; % lines; a --- header further back than this is not treated as enclosing

diary off
fid=fopen(diaryfilename,'r');
if fid<0
    diary(diaryfilename)
    fprintf('\nDiscSummary: could not reopen the diary at %s, so no summary was produced \n',diaryfilename)
    return
end
txt=textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
lines=txt{1};
diary(diaryfilename)

% STOP AT A SUMMARY THAT IS ALREADY THERE. In the normal flow this does nothing, because the summary
% is appended after the diary is read. It matters when this is called a second time on a finished
% diary, which is an easy thing to do while working on the bank: the summary echoes each failing
% line verbatim, so a second pass counts every one of them again and reports twice the failures.
sumstart=find(~cellfun(@isempty,strfind(lines,'SUMMARY OF THE WHOLE RUN')),1,'first');
if ~isempty(sumstart)
    lines=lines(1:sumstart-1);
end

nchk=0; nfail=0;
block='(before P0)'; banner='(before P0)'; section=''; sectionline=-Inf;
faillines=cell(0,1); failblock=cell(0,1); failsection=cell(0,1); failno=[];

for l_c=1:length(lines)
    l=lines{l_c};
    tk=regexp(l,'=====+\s*(P\d)[:\s]','tokens','once');
    if ~isempty(tk)
        block=tk{1};
        banner=strtrim(regexprep(strtrim(l),'^=+\s*|\s*=+$','')); % the banner without its rule of equals signs
        section=''; sectionline=-Inf; % a new subcode has started, so the previous --- header no longer encloses anything
    end
    if length(l)>10 && strncmp(strtrim(l),'---',3)
        section=strtrim(l); sectionline=l_c;
    end
    bad=0; isc=0;

    % form 1: a flag that should be one
    t=regexp(l,'this should be one:\s*(-?[\d.eE+-]+)\s*$','tokens','once');
    if ~isempty(t)
        isc=1; bad=(str2double(t{1})==0);
    end
    % form 2: should be zero - a flag if written as a bare integer, otherwise a residual
    if isc==0
        t=regexp(l,'this should be zero:\s*(-?[\d.eE+-]+)\s*$','tokens','once');
        if ~isempty(t)
            isc=1; v=abs(str2double(t{1}));
            if isempty(regexp(t{1},'[.eE]','once'))
                bad=(v~=0);
            else
                bad=(v>restol);
            end
        end
    end
    % form 3: an explicit threshold, so no convention is needed
    if isc==0
        t=regexp(l,'this should be below\s*([\d.eE+-]+):\s*(-?[\d.eE+-]+)\s*$','tokens','once');
        if ~isempty(t)
            isc=1; bad=(abs(str2double(t{2}))>str2double(t{1}));
        end
    end
    % form 4: the Monte Carlo checks, which use their own syntax
    if isc==0
        t=regexp(l,'\[MC[^\]]*:\s*([01])\]','tokens','once');
        if ~isempty(t)
            isc=1; bad=(str2double(t{1})==0);
        end
    end
    % form 5: lines that are only ever printed on the bad branch
    if isc==0 && ~isempty(strfind(l,'this should not happen'))
        isc=1; bad=1;
    end

    if isc==1
        nchk=nchk+1;
        if bad
            nfail=nfail+1;
            faillines{end+1,1}=strtrim(l); %#ok<AGROW>
            failblock{end+1,1}=block; %#ok<AGROW>
            if ~isempty(section) && (l_c-sectionline)<=seclookback
                failsection{end+1,1}=[banner,'  /  ',section]; %#ok<AGROW>
            else
                failsection{end+1,1}=banner; %#ok<AGROW>
            end
            failno(end+1,1)=l_c; %#ok<AGROW>
        end
    end
end

fprintf('\n\n')
fprintf('=================================================================================== \n')
fprintf('SUMMARY OF THE WHOLE RUN \n')
fprintf('=================================================================================== \n')
fprintf('%i checks, %i failed. Residual bar for "should be zero" lines written in e-notation: %g \n',nchk,nfail,restol)
if nfail==0
    fprintf('ALL CHECKS PASSED. \n')
else
    fprintf('\nThe %i failures, with the diary line number and the section they sit in: \n',nfail)
    for f_c=1:nfail
        fprintf('\n  [%s] diary line %i \n',failblock{f_c},failno(f_c))
        if ~isempty(failsection{f_c})
            fprintf('  in section: %s \n',failsection{f_c})
        end
        fprintf('  %s \n',faillines{f_c})
    end
    fprintf('\nA failure here is not automatically a toolkit bug. A [T0] or [T1] failure is an identity \n')
    fprintf('that did not hold and is worth chasing. A [T2] failure is a measured accuracy number \n')
    fprintf('against a bar this bank chose, so it can equally mean the bar is wrong, or that the \n')
    fprintf('method simply cannot do better - read the section the line sits in before the number. \n')
end
fprintf('=================================================================================== \n')

end
