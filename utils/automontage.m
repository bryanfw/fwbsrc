function automontage(img, range )

% pass real values 3d image, where the image information is in dims 1 & 2
% optionally pass range as a [min max] array

if nargin < 2
    range(1) = min(0,min(img(:))); 
    range(2) = max(img(:)); 
end

montage(permute(img,[1 2 4 3]), [range(1), range(2)])