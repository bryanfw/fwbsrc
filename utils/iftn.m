function y = iftn(x)
% Performs ifftshift(ifftn(fftshift(input)))
% 2D inverse FT

y = ifftshift(ifftn(fftshift(x))); 

