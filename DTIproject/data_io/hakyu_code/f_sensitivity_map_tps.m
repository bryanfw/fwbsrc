
function [Sout,BWout] = f_sensitivity_map_tps(I_body_m, I_ref_3d, vox_dilate)
% [f_sensitivity_map_tps] generates sensitivity map using 2D
% thin-plate spline.
%
% INPUT
%   I_body_m:   Bodycoil image, for one slice
%   I_ref_3d:   Individual coil image, for one slice for all coils
%   vox_dilate: Number of voxels to dilate to generated rim region
%
% OUTPUT
%   S:  Coil sensitivity map, for one slice and for all coils
%   BW: Mask used for generating the coil sensitivity map
%
%
% Last modified
% 2010.09.28.
%   Generated from [test_f_sensitivity_map_tps.m].
% 2010.09.30.
%   Use BW dilation for generating final sensitivity map and mask.
% 2012.05.16.
%   Remove text output for showing computation time.
% 
% Ha-Kyu



%% Preliminary
[Ny,Nx,Nc] = size(I_ref_3d); %[M,P,coil] direction



%% Generate BW
I = abs(I_body_m);
I = imadjust(I);
BW0 = im2bw(I,graythresh(I)*0.5);
BW0 = imfill(BW0,'holes');

% Find max area which supposed to be the main object.
L = bwlabel(BW0,8);
stats = regionprops(L,'Area');
idx = find([stats.Area]==max([stats.Area]));
BW = ismember(bwlabel(BW0,8),idx);
%[ind_row,ind_col] = find(BW);
%min_ind_col = min(ind_col);

% Dilate BW.
Rfit = vox_dilate;             % default:3
se = strel('disk',Rfit,0);
BWdil = imdilate(BW,se);

% Clear.
clear  BW0  L



%% Generate random sampling points and mask

% Random sampling points in the original image space.
nrand = round(Ny*Nx*0.1); % 10% of total voxels
xrand_v = rand(nrand,1); xrand_v = round((Nx-1)*xrand_v+1);
yrand_v = rand(nrand,1); yrand_v = round((Ny-1)*yrand_v+1);

% Generate mask with the random sampling points.
ind_rand_v = sub2ind([Ny,Nx],yrand_v,xrand_v);
BW_rand_m = zeros(Ny,Nx);
BW_rand_m(ind_rand_v) = 1;
BW_rand_mask_m = BW_rand_m.*BW;

% Find min and max row and col sampling points.
[ind_rand_row_v,ind_rand_col_v] = find(BW_rand_mask_m);
ind_min_row = min(ind_rand_row_v);
ind_max_row = max(ind_rand_row_v);
ind_min_col = min(ind_rand_col_v);
ind_max_col = max(ind_rand_col_v);

% Find fit region for TPS.
% I_body_m and I_ref_3d always have M_ORI along row (Ny) direction.
ind_y_start = ind_min_row-vox_dilate-5; % to cover BWdil, subtract(add) 5 voxels
ind_y_end = ind_max_row+vox_dilate+5;
ind_x_start = ind_min_col-ind_min_col+1;
ind_x_end = ind_max_col-ind_max_col+Nx;
if ind_y_start < 1
    ind_y_start = 1;
end
if ind_y_end > Ny
    ind_y_end = Ny;
end
Nxr = length(ind_x_start:ind_x_end); % same as Nx
Nyr = length(ind_y_start:ind_y_end); % reduced from Ny

% Generate mask for TPS fit region.
% When TPS fit region is too big, e.g., [Ny,Nx], the coil sensitivity map
% can be wrong in which the sensitivity profile doesn't match the raw
% sensitivity map, here raw sensitivity map is ref_coil / bodycoil.
BW_rand_fit_m = BW_rand_mask_m(ind_y_start:ind_y_end,ind_x_start:ind_x_end);

% Generate reduced mask including BW_rand_fit_m.
% BWr = BW(ind_y_start:ind_y_end,ind_x_start:ind_x_end);

% Get random sampling coordinates inside BW_rand_fit_m (measured points).
[rand_row_idx_v,rand_col_idx_v] = find(BW_rand_fit_m);
rand_idx_v = find(BW_rand_fit_m);
xy = [rand_col_idx_v'; rand_row_idx_v'];
clear  BW_rand_m  xrand_v  yrand_v  ind_rand_v  BW_rand_m



%% Generate raw coil sensitivity map
A = I_body_m;
smap_3d = I_ref_3d./repmat(A,[1,1,Nc]);
smap_thresh_3d = smap_3d.*repmat(BW,[1 1 Nc]);
clear  smap_3d  A



%% Generate interpolation and grid points and M, K, P matrices

% Grid points.
[x_m,y_m] = meshgrid(1:Nxr,1:Nyr);
% [x_m,y_m] = meshgrid(1:Nx,1:Ny); % too wide fit region spoils coil sensitivity
xyg = [x_m(:)'; y_m(:)'];
clear  x_m  y_m


%--------------------------------------------------------
% Calculate K, P, O for interpolation-measurement points.
% This part is to calculate coefficients.

% Generate interpolation points.
% x_v = min(xy(1,:)):max(xy(1,:));
% y_v = min(xy(2,:)):max(xy(2,:));
% [x_m,y_m] = meshgrid(x_v,y_v);
% xyi = [x_m(:)'; y_m(:)'];
xyi = xy; % #measurement points == #interpolation points
clear  x_m  y_m  x_v  y_v

% Calculate K.
a = xy(1,:)'+1i*xy(2,:)'; % measured
b = xyi(1,:)'+1i*xyi(2,:)'; % interp
a1 = repmat(a,1,length(b));
b1 = repmat(b.',length(a),1);
clear  a  b

d = a1-b1;  clear  a1  b1
R = d.*conj(d);  clear  d
ind1 = (R==0);
ind2 = (R>0);
K = zeros(size(R));
K(ind1) = 0;
K(ind2) = log(R(ind2)).*R(ind2);
clear  R

% Calculate P.
P = [ones(1,length(xy)); xy];

% Calculate O.
O = zeros(3,size(K,2)+3-size(K,1));

% Calculate M.
M = [K, P.'; P, O];
clear  K  P  O


%----------------------------------------------
% Calculate K, P for grid-interpolation points.
% This part is to calculate TPS function.

% Calculate K.
a = xyg(1,:)'+1i*xyg(2,:)';
b = xyi(1,:)'+1i*xyi(2,:)';
a1 = repmat(a,1,length(b));
b1 = repmat(b.',length(a),1);
clear  a  b

d = a1-b1;  clear  a1  b1
R = d.*conj(d);  clear  d
ind1 = (R==0);
ind2 = (R>0);
K = zeros(size(R));
K(ind1) = 0;
K(ind2) = log(R(ind2)).*R(ind2);
clear  R

% Calculate P.
P = [ones(1,size(xyg,2)); xyg];



%% Fit using TPS

% Reserve output.
smap_fit_3d = zeros(Ny,Nx,Nc,'single');

% TPS fit.
for ind_coil = 1:Nc
    tic    
    
    % Fit SMAP for real and imaginary values.    
    for ind_typ = 1:2 % real and imaginary
        if ind_typ==1
            smap_m = real(smap_thresh_3d(ind_y_start:ind_y_end, ...
                ind_x_start:ind_x_end,ind_coil));
        else
            smap_m = imag(smap_thresh_3d(ind_y_start:ind_y_end, ...
                ind_x_start:ind_x_end,ind_coil));
        end
        smap_val_v = smap_m(rand_idx_v);
        clear  smap_m
        
        % Y.
        y = [smap_val_v; zeros(3,1)];
        
        % Calculate alpha and beta.
        [Q,R] = qr(M);
        y = Q'*y;
        coeff = R \ y;
        %coeff = M \ y;
        clear  Q  R  y % keep M
        
        
        %--------------------------------
        % Calculate TPS function.
                
        % Calculate TPS.
        f = [K,P']*coeff;
        z_m = reshape(f,Nyr,Nxr);
        %z_m = reshape(f,Ny,Nx);    % too wide fit region messes coil sensitivity
        clear  f % keep K, P
        
        % Reserve fitted SMAP.
        if ind_typ==1
            smap_fit_real_m = z_m;
        else
            smap_fit_imag_m = z_m;
        end
        clear  z_m
    end
    smap_fit_m = complex(smap_fit_real_m,smap_fit_imag_m);
    
    t=toc;
    %fprintf('    SMAP for coil[%d] calculated in [%.3f]sec\n',ind_coil,t)
    
    % Output.
    smap_fit_pad_m = zeros(Ny,Nx);
    smap_fit_pad_m(ind_y_start:ind_y_end,ind_x_start:ind_x_end) = smap_fit_m;
    %smap_fit_3d(:,:,ind_coil) = smap_fit_pad_m;
    smap_fit_3d(:,:,ind_coil) = smap_fit_pad_m.*BWdil;
    
    %smap_fit_3d(ind_y_start:ind_y_end,ind_x_start:ind_x_end,ind_coil) = smap_fit_pad_m;
    %smap_fit_3d(:,:,ind_coil) = smap_fit_m; % too wide fit region messes coil sensitivity
    
    % Clear data.
    clear  smap_fit_real_m  smap_fit_imag_m  smap_fit_m  smap_fit_pad_m
end



%% Output
Sout = smap_fit_3d;
BWout = BWdil;

% Clear and pack.
clear  smap_fit_3d  I_ref_3d  I_body_m  BWdil  BW_*  
pack

%fprintf('        S generated.\n')


%% END






