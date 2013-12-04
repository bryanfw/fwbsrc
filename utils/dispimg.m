function dispimg(im,range)

if nargin < 2
    range(1) = min(0,min(im(:))); 
    range(2) = max(im(:)); 
end

imshow(im, [range(1), range(2)]); 

end