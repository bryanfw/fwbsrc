function y = ift2(x)

% Performs ifftshift(ifft2(fftshift(input)))
% 2D inverse FT

if size(x,3) ==1;
    y = ifftshift(ifft2(ifftshift(x))); 
else
    y = zeros(size(x));
    for i=1:size(x,3);
        y(:,:,i) = ift2(x(:,:,i));
    end
end
