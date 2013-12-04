function y = ft(x,dim)
if nargin<2
    dim = 1;
end
% Performs fftshift(fft(ifftshift(input)))
% forward FT
y = fftshift(fft(ifftshift(x,dim),[],dim),dim); 

