function y = ft2(x)
% Performs fftshift(fft2(fftshift(input)))
% 2D forward FT

if ismatrix(x);
    y = fftshift(fft2(ifftshift(x))); 
else
    y = zeros(size(x));
    for i=1:size(x,3);
        y(:,:,i) = ft2(x(:,:,i));
    end
end

