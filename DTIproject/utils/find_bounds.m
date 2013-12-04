function [mn,mx] = find_bounds(image)
% image = (image-min(image(:)))/(max(image(:))-min(image(:)));
mn = prctile(image(:),1);
mx = prctile(image(:),99);