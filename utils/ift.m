function y = ift(x,dim)
% Performs ifftshift(ifft(fftshift(input)))
% inverse FT
if nargin<2
    dim = 1;
end
y = ifftshift(ifft(fftshift(x,dim),[],dim),dim); 

