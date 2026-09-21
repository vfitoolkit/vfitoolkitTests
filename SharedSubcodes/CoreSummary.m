function CoreSummary(diaryfilename)
% Read the diary of a Core* test bank back and print one verdict for the whole run.
%
% This is the Core-bank counterpart of DiscSummary (DiscretizationMethodTests/DiscSummary.m).
%
% WHY THIS EXISTS. A Core bank prints on the order of a thousand checks into a diary of several
% thousand lines, and every one of them is a number that has to be compared against a bar in the
% reader's head. "Did the run pass" is not a question anyone should answer by scrolling, and the
% check COUNT is worth as much as the verdict: a check that silently stopped being executed - a
% call commented out, a tier skipped, a subcode that returned early - shows up as a smaller count
% and as nothing else at all.
%
% The diary is closed, read, and reopened in append mode, so this summary lands in the diary itself.
%
% THE THREE FORMS, and why they are not held to the same bar:
%
%   'this should be zero: 2.665e-15'          IDENTITY. Two ways of computing the same object must
%                                             agree. Held to restol.
%   'should give zero: 0.000e+00'             IDENTITY too (the QH beta0=1 and EZ riskaversion=0
%                                             degenerate-parameter checks). Same bar.
%   'V   should be ~0: 0.00000001'           IDENTITY, in the InfHorz VFI-algorithm banks' spelling.
%   'Pol should be  0: 0.00000000'           IDENTITY too (an exact Policy match). Same bar.
%                                            NB those banks print %.8f, which cannot tell zero from
%                                            a few ULP - see [[test-bank-zero-checks-ulp-floor]].
%   'this should be close to zero: 4.839e-03' NOT an identity, and the wording is the bank telling
%                                             you so. It is the agent distribution with and without
%                                             the grid interpolation layer, which genuinely differ:
%                                             interpolation puts mass between grid points. On the
%                                             ExpAssetU diary of 2026-09-19 these run to 2.5e-2
%                                             while every identity on the same run is under 1e-8.
%                                             Held to closetol, and always reported, because the
%                                             number drifting is informative even when it passes.
%
% RESTOL, AND WHERE IT COMES FROM. restol cannot be set as tight as DiscSummary's 1e-10. These
% banks compare value functions, and V reaches about 1e7 at the poor corner of the grid, so one
% unit in the last place there is about 4e-9: a difference of a few ULP is not a disagreement, it
% is two orders of operations on the same arithmetic.
%
% The bar was not chosen by feel. Across the three diaries current on 2026-09-20 (CoreFHorzTests,
% CoreFHorzExpAssetTests, CoreFHorzExpAssetUTests - 3031 identity checks in total) the residuals
% are sharply bimodal, and there is an empty decade to put the bar in:
%    exactly zero      2524
%    below 1e-12        468     ordinary double-precision noise
%    1e-12 to 1e-10       0     <- empty
%    1e-10 to 1e-8       39     the V~1e7 ULP band: ValueFnFromPolicy with grid interpolation
%    1e-8 and above       0     <- empty; the largest residual anywhere is 5.588e-09
% So restol=1e-7 sits in empty space, a factor of 18 above the largest residual any of these banks
% has ever printed and five orders below the scale of the 'close to zero' checks. Re-derive it if
% a bank starts printing residuals into that gap; do not quietly raise it to make a run pass.
%
% The passing identities are split three ways in the report - exactly zero, below noisetol, and
% between noisetol and the bar - so that drift INTO the ULP band is visible long before anything
% fails. A failing identity is counted as a failure and in none of the three bands.
%
% NAMES CONTAIN DIGITS. Values are matched only where a number stands on its own, because the
% cross-tests label their values with names like 'AllStats.a2.Mean' and 'a1_2' - reading the digits
% out of those turns a line of exact zeros into a failure.
%
% NEARTOL, AND ITS OWN EMPTY GAP. 'should be near zero' is the CoreStationaryGeneralEqm bank's
% wording for comparing two GE solves, and its values do not live on the same scale as either an
% identity or a distribution difference, so it gets its own bar. On that bank's diary the plain
% 'near zero' values are 1.1e-16 to 7.6e-5 (transform round trips, then solver-accuracy
% comparisons) and then two outliers at 6.5e-2 and 6.9e-2 - a gap of a factor of 854. neartol=1e-3
% sits in it, 13x above the highest ordinary value and 65x below the outliers. Those two outliers
% are 'fminalgo 1 vs 8', i.e. two optimizers disagreeing about the same equilibrium by about 7%,
% and they SHOULD be flagged: a bar loose enough to pass them would hide the thing worth seeing.
% Note the bank marks its own known-loose checks '(loose: <reason>)', and those are taken at
% closetol instead - it is saying so itself.
%
% A LINE MAY CARRY SEVERAL CHECKS. The cross-tests print
%    'Cross test (noa1): this should be zero: V 0.000e+00, Policy 0.000e+00, Dist 6.939e-18'
% which is three checks, and is counted as three. Anything in parentheses after the value is
% diagnostic detail, not a further check ('(rel 1.3e-16, max|V|=1.1e+07, worst at [3] of [5])',
% '(2 of 100 entries differ)'), so the line is cut at the first '(' after the colon.
%
% READING A POLICY FAILURE. A 'Divide-and-conquer (DC2A)' or similar line reporting exactly 1 or 2,
% on a bank whose V agrees to the ULP on the neighbouring line, is usually not a broken identity: it
% is the known DC2A tie artefact. Where two decisions give exactly the same payoff, the narrow-band
% divide-and-conquer search breaks the tie differently from the full search, so the Policy INDEX
% differs by one grid point while the value function does not differ at all. It was confirmed benign
% on the GPU for fig 42 of CoreFHorzExpAssetTests on 2026-09-19 (6 exact ties, d2 payoff-irrelevant)
% and the bank fix was declined. These are still reported rather than passed, because the count and
% location moving IS information; they are simply not news.
%
% PER PART. Failures and counts are attributed to the doPart block they sit in, using the
% '===== doPart(k): label =====' banner each part prints. A bank with no such banners still works;
% everything is then attributed to '(whole run)'.

restol=1e-7;    % identity bar, set in the empty decade above the largest observed residual (see above)
closetol=0.1;   % 'close to zero' bar; largest seen across the three banks as of 2026-09-20 is 2.5e-2
noisetol=1e-10; % splits ordinary double noise from the V~1e7 ULP band; reporting only, not a bar
neartol=1e-3;   % 'should be near zero' bar, set in that population's own empty decade (see above)

diary off
fid=fopen(diaryfilename,'r');
if fid<0
    diary(diaryfilename)
    fprintf('\nCoreSummary: could not reopen the diary at %s, so no summary was produced \n',diaryfilename)
    return
end
txt=textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
lines=txt{1};
diary(diaryfilename)

% STOP AT A SUMMARY THAT IS ALREADY THERE, so that calling this a second time on a finished diary
% (an easy thing to do while working on a bank) does not count the echoed failure lines again.
sumstart=find(contains(lines,'SUMMARY OF THE WHOLE RUN'),1,'first');
if ~isempty(sumstart)
    lines=lines(1:sumstart-1);
end

part='(whole run)';
partnames=cell(0,1); partchk=[]; partfail=[];
nchk=0; nfail=0; nexact=0; nnoise=0; nulp=0; nloose=0; nother=0; nnear=0;
maxnear=0; maxnearline=''; maxnearno=0;
maxres=0; maxresline=''; maxresno=0;
maxloose=0; maxlooseline=''; maxlooseno=0;
faillines=cell(0,1); failpart=cell(0,1); failno=[];

for l_c=1:length(lines)
    l=lines{l_c};

    tk=regexp(l,'=====+\s*doPart\((\d+)\)\s*:\s*(.*?)\s*=====','tokens','once');
    if ~isempty(tk)
        part=['doPart(',tk{1},'): ',tk{2}];
    end

    % Which form is this line, if any? 'close to zero' must be tested BEFORE 'be zero', since the
    % looser text contains the stricter pattern's words.
    % WHICH FORM IS THIS LINE? Order matters: the looser wordings contain the stricter ones.
    isloose=0; isflag=0; isinv=0; isnear=0;
    % Loose: 'should be close to zero:', 'should be very similar (small):', 'should be small:', and
    % 'should be near zero (loose: <reason>)' - that last one is a bank marking its OWN check as
    % deliberately loose, which is worth honouring, and it must be tested before the plain
    % 'near zero' below or it would be held to the tighter bar. The InfHorz TPath banks write the
    % first inside a parenthesis - '(these should be close to zero): Capital=...' - so text is
    % allowed between the wording and the colon; without that the line is invisible.
    rest=regexp(l,'should be (?:close to zero|very similar|small|near zero \(loose)[^:]*:\s*(.*)$','tokens','once');
    if ~isempty(rest)
        isloose=1;
    end
    % Inverted: 'should be WELL ABOVE zero:' and 'should NOT be zero:' are passed by being far FROM
    % zero, so they are checked the other way up. Tested before the identity form, which would
    % otherwise read them backwards.
    if isempty(rest)
        rest=regexp(l,'should (?:be WELL ABOVE zero|NOT be zero)[^:]*:\s*(.*)$','tokens','once');
        if ~isempty(rest)
            isinv=1;
        end
    end
    % Near: the plain 'should be near zero:'. Its own bar - see neartol above.
    if isempty(rest)
        rest=regexp(l,'should be near zero[^:]*:\s*(.*)$','tokens','once');
        if ~isempty(rest)
            isnear=1;
        end
    end
    % Flag: 'this should be one: 1' fails at 0. DiscSummary carries the same form.
    if isempty(rest)
        rest=regexp(l,'should be one\s*:\s*(.*)$','tokens','once');
        if ~isempty(rest)
            isflag=1;
        end
    end
    % Identity: 'should be zero:', 'should give zero:', and two spellings of the same claim in the
    % InfHorz banks - 'should be ~0:' (a float residual) and 'should be  0:' (an exact Policy
    % match). Text is allowed between the wording and the colon, because the TPath banks label the
    % value there: 'this should be zero, Policy: 6.000e+00'. Without that the Policy comparison on a
    % constant transition path was invisible, and CoreInfHorzVFIAlgoTests reported 0 checks in 1929
    % lines of them.
    if isempty(rest)
        rest=regexp(l,'should (?:be|give) (?:zero|~?\s*0)[^:]*:\s*(.*)$','tokens','once');
    end
    if isempty(rest)
        continue
    end

    % WHERE THE VALUE IS. Normally it is everything after the colon that follows the wording. Two
    % different things break that, and they have to be told apart:
    %   '... should be zero, AgentDist (note: tolerance=1e-6): 4.907e-06'   ONE value; the colon we
    %       split on sits INSIDE an annotation, so the '1e-6' is not a check.
    %   '... should be small, r: 0.0727387574, w: 0.3281381734'             TWO values, both checks.
    % The tell is an unbalanced ')' in what follows: it means the split colon was inside a
    % parenthetical that opened before it, so the value is after the LAST colon on the line.
    % Otherwise the text is a list of labelled values and all of them count. Getting this wrong
    % either invents a check out of a quoted tolerance, or silently drops the first of two.
    s=rest{1};
    if sum(s==')')>sum(s=='(')
        ic=find(l==':',1,'last');
        s=l(ic+1:end);
    end
    % Then cut any diagnostic parenthetical that trails the value.
    ip=strfind(s,'(');
    if ~isempty(ip)
        s=s(1:ip(1)-1);
    end
    % A value is a number standing on its own, never digits embedded in a name. Without the
    % lookbehind, 'AllStats.a2.Mean 0.000e+00' yields the 2 of 'a2' as a value of 2, and every one
    % of the RiskyAsset ExpAssetu cross-tests reads as a failure while printing four exact zeros.
    v=str2double(regexp(s,'(?<![A-Za-z_0-9.])[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?','match'));
    v=v(~isnan(v));
    if isempty(v)
        continue
    end
    v=abs(v);

    % Record this part if it is the first check seen under it
    ip=find(strcmp(partnames,part),1);
    if isempty(ip)
        partnames{end+1,1}=part; partchk(end+1,1)=0; partfail(end+1,1)=0; %#ok<AGROW>
        ip=numel(partnames);
    end

    bar=restol;
    if isloose, bar=closetol; elseif isnear, bar=neartol; end
    for v_c=1:numel(v)
        nchk=nchk+1; partchk(ip)=partchk(ip)+1; %#ok<AGROW>
        if isflag || isinv
            nother=nother+1;
        elseif isnear
            nnear=nnear+1;
            if v(v_c)>maxnear
                maxnear=v(v_c); maxnearline=strtrim(l); maxnearno=l_c;
            end
        elseif isloose
            nloose=nloose+1;
            if v(v_c)>maxloose
                maxloose=v(v_c); maxlooseline=strtrim(l); maxlooseno=l_c;
            end
        else
            % The three bands partition the identities that PASS; a failure is counted as a
            % failure and nowhere else, so the bands cannot be padded out by bad numbers.
            if v(v_c)==0
                nexact=nexact+1;
            elseif v(v_c)<noisetol
                nnoise=nnoise+1;
            elseif v(v_c)<=restol
                nulp=nulp+1;
            end
            if v(v_c)>maxres
                maxres=v(v_c); maxresline=strtrim(l); maxresno=l_c;
            end
        end
        if isflag
            isbad=(v(v_c)==0);            % a flag fails at zero
        elseif isinv
            isbad=(v(v_c)<=restol);       % 'WELL ABOVE zero' fails by being near zero
        else
            isbad=(v(v_c)>bar);
        end
        if isbad
            nfail=nfail+1; partfail(ip)=partfail(ip)+1; %#ok<AGROW>
            faillines{end+1,1}=strtrim(l); %#ok<AGROW>
            failpart{end+1,1}=part; %#ok<AGROW>
            failno(end+1,1)=l_c; %#ok<AGROW>
        end
    end
end

fprintf('\n\n')
fprintf('=================================================================================== \n')
fprintf('SUMMARY OF THE WHOLE RUN \n')
fprintf('=================================================================================== \n')
fprintf('%i checks, %i failed. \n',nchk,nfail)
fprintf('  identities ("should be zero", bar %g): %i \n',restol,nchk-nloose-nother-nnear)
fprintf('    passing: %i exactly zero, %i below %g (double noise), %i in the V~1e7 ULP band above it \n',nexact,nnoise,noisetol,nulp)
if nchk-nloose-nother-nnear>0 && maxresno>0
    fprintf('    largest: %.3e  (diary line %i) \n',maxres,maxresno)
    fprintf('      %s \n',maxresline)
elseif nchk-nloose-nother-nnear>0
    fprintf('    every identity is an exact zero \n')
end
fprintf('  "close to zero" checks (bar %g, NOT identities - see the header of CoreSummary.m): %i \n',closetol,nloose)
if nloose>0
    fprintf('    largest: %.3e  (diary line %i) \n',maxloose,maxlooseno)
    fprintf('      %s \n',maxlooseline)
end
if nnear>0
    fprintf('  "near zero" checks (bar %g): %i \n',neartol,nnear)
    fprintf('    largest: %.3e  (diary line %i) \n',maxnear,maxnearno)
    fprintf('      %s \n',maxnearline)
end
if nother>0
    fprintf('  flags ("should be one") and inverted ("should NOT be zero") checks: %i \n',nother)
end

if ~isempty(partnames)
    fprintf('\nChecks by part (only the parts that actually ran appear here): \n')
    for p_c=1:numel(partnames)
        fprintf('  %-52s %5i checks, %i failed \n',partnames{p_c},partchk(p_c),partfail(p_c))
    end
end

if nfail==0
    fprintf('\nALL CHECKS PASSED. \n')
else
    fprintf('\nThe %i failures, with the diary line number and the part they sit in: \n',nfail)
    for f_c=1:nfail
        fprintf('\n  [%s] diary line %i \n',failpart{f_c},failno(f_c))
        fprintf('  %s \n',faillines{f_c})
    end
    fprintf('\nBefore chasing one: check the COUNT above against the last run of this bank. A failure \n')
    fprintf('is only meaningful if the same checks ran. And a value just over %g may be the ULP floor \n',restol)
    fprintf('moving rather than an identity breaking - the largest-residual line above says how much \n')
    fprintf('headroom the bar had. \n')
end
fprintf('=================================================================================== \n')

end
