function y = ftn(x)

% Performs fftshift(fftn(fftshift(input)))
% 2D forward FT

y = fftshift(fftn(ifftshift(x))); 

