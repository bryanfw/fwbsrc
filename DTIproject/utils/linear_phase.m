function ph = linear_phase(sz, shift)
% LINEAR_PHASE creates a linear phase (complex expontial of magnitude 1) that 
% will shift your image the input number of pixels
%
% phase = LINEAR_PHASE(size_vector, shift_vector)
% size_vector is two element array of the same type as given by size()
% shift_vector is a two element array giving +/- the number of voxel shift
%   desired. Negative values cause upward leftward movement.
%
% Only designed with 2D in mind.
%

% Frederick Bryan, Vanderbilt, 2013

shiftx = shift(1);
shifty = shift(2);

[phx,phy] = meshgrid((-floor(sz(2)/2):floor(sz(2)/2)-1)/sz(2),...
    (-floor(sz(1)/2):1:floor(sz(1)/2)-1)/sz(1));
ph = 2*pi*(phx*-shiftx+phy*shifty); 

% 
% 
% rows = sz(1);
% cols = sz(2);
% 
% [phy,phx] = ndgrid((0:floor(rows-1))/rows, ...
%     (-floor(cols/2):floor(cols/2-1))/cols ); 
% 
% ph = 2*pi*(phy)*pix_shift;


