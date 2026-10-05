function [paths,values]=CoreFHorzPType_StructLeaves(S,prefix,paths,values)
% List every numeric/logical leaf of a (nested) structure: paths{k} is the dotted field path (e.g. '.assets.low.Mean'), values{k}
% the gathered array. Cell arrays and character fields (e.g. .Notes) are skipped. Used by the whichcombos tests to compare a
% selective run against a full run leaf by leaf.

if nargin<2
    prefix='';
    paths={};
    values={};
end
fnames=fieldnames(S);
for k=1:numel(fnames)
    x=S.(fnames{k});
    p=[prefix,'.',fnames{k}];
    if isstruct(x)
        if isscalar(x)
            [paths,values]=CoreFHorzPType_StructLeaves(x,p,paths,values);
        end
    elseif isnumeric(x) || islogical(x)
        paths{end+1}=p; %#ok<AGROW>
        values{end+1}=gather(x); %#ok<AGROW>
    end
end

end
