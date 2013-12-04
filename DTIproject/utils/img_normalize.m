function [out] =img_normalize(in,prctiles,mask)
% out = IMG_NORMALIZE(in)
% utility function used to normalize intensity to between 0.5 and 99.5
% percentile and remap everyting [0 1]

out = in;

if nargin < 2 || isempty(prctiles)
    prctiles = [.5 99.5];
end
if nargin==3 % if mask given, mask image
    if isa(mask,'logical') 
        in = in(mask);
    elseif all(ismember(mask(:),[0 1])); % check if all elements are in {0,1}
        mask = logical(mask);
        in = in(mask);
    else
        error('mask elements must be [0 1] if not class:logical');
    end
end

% cast to double
in = double(in);

% saturate btwen .5 and 99.5 percentile
lowval = prctile(in(:), prctiles(1));
highval = prctile(in(:), prctiles(2));

out(out<lowval) = lowval;
out(out>highval) = highval;

% scale to one;
out = (out-lowval)/(highval-lowval);

