%
% SCROLLY - scroll a 2D data matrix in the Y (second dimension) by a known number of pixels
%
% scroll_data = scrollY(data,pixels)
%
% data                 - input 2D data matrix 
% pixels               - number of pixels to scroll in Y
% scroll_data          - output 2D data matrix
%
%

function sp = scrollY(data,pixels);

data_size = size(data);
x_size    = data_size(1);
y_size    = data_size(2);

for I = 1:y_size
    ind     = 1 + rem(I + y_size - pixels + 1, y_size);
    sp(:,I) = data(:,ind);
end