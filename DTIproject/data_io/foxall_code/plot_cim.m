%
% PLOT_CIM  - Show the corrected IEPI image
%


figure

imagesc(abs(cpc_image));
colormap(gray(256))
axis('image')
