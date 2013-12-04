%
% SCROLLX - scroll a 2D data matrix in the X (first dimension) by a known number of pixels
%
% scroll_data = scrollX(data,pixels)
%
% data                 - input 2D data matrix 
% pixels               - number of pixels to scroll in X
% scroll_data          - output 2D data matrix
%
%

function sp = scrollX(data,pixels);

data_size = size(data);
x_size    = data_size(1);
y_size    = data_size(2);

for I = 1:x_size
    ind     = 1 + rem(I + x_size - pixels + 1, x_size);
    sp(I,:) = data(ind,:);
end