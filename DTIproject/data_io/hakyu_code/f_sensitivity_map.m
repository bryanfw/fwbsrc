
function [S,BW] = f_sensitivity_map(I_body_m, I_ref_3d, fitorder, fitsize_out, fitsize_in, vox_dilate)
% [f_sensitivity_map] generates sensitivity map.
%
% INPUT
%   I_body_m:   Bodycoil image, for one slice
%   I_ref_3d:   Individual coil image, for one slice for all coils
%   fitorder:   2D polynomial fitting order
%   fitsize_out:    Size of support region for fitting outside of rim region
%   fitsize_in:     Size of support region for fitting inside of rim region
%   vox_dilate:     Number of voxels to dilate to generated rim region
%
% OUTPUT
%   S:  Coil sensitivity map, for one slice for all coils
%   BW: Mask used for generating the coil sensitivity map
%
% NOTE
%   The 'rim' region is the region dilated radially from the original
%   object mask. See the cell '%% Generate BW'.
%
%
% Last modified
% 2009.07.30. v4. Generate coil sensitivity map using Pruessmann's original
%                 method (MRM 1999). This is tested by
%                 [test_gen_smap_v2.m].
%               
%                   fitorder,       default:2
%                   fitsize_out,    default:10 (15)
%                   fitsize_in,     default:5 (9)
%                   vox_dilate,     default:3 (6)
%
% 2009.09.26.     v5 is generated using smoothing (whole) and 2D fitting
%                 (rim).
%                 v5 doesn't look better than this. Just use this version
%                 v4.
% 2009.10.21.     Change default value for fitsize_out (10->15), fitsize_in
%                 (5->9), vox_dilate(3->6) to cover distorted object
%                 regions due to field inhomogeneity.
% 2009.11.18.     Generated from [sensitivity_map_v4.m]. This doesn't
%                 generate PSI because PSI is calculated in
%                 [REF_readsaveshow_scan********.m] using NOI data.
% 2010.08.09.
%   This is generated from [sensitivity_map_v6.m].
%   This process can be accelerated in the future, if several voxels are
%   fitted, not for single voxel as used here.



%% Preliminary
[Ny,Nx,Nc] = size(I_ref_3d);


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

% Dilate BW.
Rfit = vox_dilate;             % default:3
se = strel('disk',Rfit,0);
BWfit = imdilate(BW,se);

% Outer rim region.
BWrim = BWfit-BW;

clear L



%% Generate smap
A = I_body_m;
Smap_3d = I_ref_3d./repmat(A,[1,1,Nc]);
Smap_thresh_3d = Smap_3d.*repmat(BW,[1 1 Nc]);



%% Generate smoothing kernel
% This procedure skipped.

%% Smooth thresholded smap
% This procedure skipped.

%% Generate low-pass filter
% This procedure skipped.

%% Low-pass filtering
% This procedure skipped.



%% 2D polynomial fit

% Fit region coordinates and indices.
[y_v,x_v] = find(BWfit);
ind_fit_v = find(BWfit);
ind_rim_v = find(BWrim);

% Fit parameters.
maxOrder = fitorder;    % default:2
n_out = fitsize_out;    % default:10, 15
n_in = fitsize_in;      % default:5, 9
coil_v = 1:Nc;

% Reserve output.
Smap_fit_3d = zeros(size(Smap_thresh_3d),'single');

% Generate coil sensitivity map.
for coil = coil_v
    tic

    Smap_m = Smap_thresh_3d(:,:,coil);
    Smap_fit_m = Smap_m*0;
    
    h = waitbar(0,'Please wait...');
    for ind_vox = 1:length(y_v)
        
        % Waitbar.
        if mod(ind_vox,100)==0
            waitbar(ind_vox/length(y_v),h,sprintf('COIL[%d], voxel[%d] of [%d]',coil,ind_vox,length(y_v)))
        end
        
        % Size of fit region.
        if any(ind_fit_v(ind_vox)==ind_rim_v) % outer rim region
            n = n_out;
        else
            n = n_in;
        end
        
        % Voxel coordinate for fit in whole image.
        r = y_v(ind_vox);
        c = x_v(ind_vox);
        
        % Define fit region coordinates.
        r_v = (r-n:r+n)';
        c_v = (c-n:c+n)';
        
        % Check out of boundary coordinates.
        indr_range_v = (r_v >= 1) & (r_v <= Ny);
        indc_range_v = (c_v >= 1) & (c_v <= Nx);
        r_v = r_v(indr_range_v);
        c_v = c_v(indc_range_v);
        
        if length(find(indr_range_v)) < length(find(indc_range_v))
            total_n = length(r_v);
            current_n = length(c_v);
            indc = find(c==c_v);            
            vox_remain_left = indc-1;
            vox_remain_right = current_n-indc;
            vox_keep_left = floor(total_n/2);
            vox_keep_right = ceil(total_n/2)-1;
            
            delta_left = vox_remain_left-vox_keep_left;
            delta_right = vox_remain_right-vox_keep_right;
            
            if delta_left <= 0
                ind_start = 1;
                ind_end = indc+vox_keep_right-delta_left;
            elseif delta_right <= 0
                ind_start = indc-vox_keep_left+delta_right;
                ind_end = current_n;
            else
                ind_start = indc-vox_keep_left;
                ind_end = indc+vox_keep_right;
            end
            if length(ind_start:ind_end)~=total_n
                error('check ind_start and ind_end of c_v.')
            end            
            c_v = c_v(ind_start:ind_end);
        else
            total_n = length(c_v);
            current_n = length(r_v);
            indr = find(r==r_v);            
            vox_remain_left = indr-1;
            vox_remain_right = current_n-indr;
            vox_keep_left = floor(total_n/2);
            vox_keep_right = ceil(total_n/2)-1;
            
            delta_left = vox_remain_left-vox_keep_left;
            delta_right = vox_remain_right-vox_keep_right;
            
            if delta_left <= 0
                ind_start = 1;
                ind_end = indr+vox_keep_right-delta_left;
            elseif delta_right <= 0
                ind_start = indr-vox_keep_left+delta_right;
                ind_end = current_n;
            else
                ind_start = indr-vox_keep_left;
                ind_end = indr+vox_keep_right;
            end
            if length(ind_start:ind_end)~=total_n
                error('check ind_start and ind_end of r_v.')
            end            
            r_v = r_v(ind_start:ind_end);
        end
        
        row_fit = find(r==r_v);
        col_fit = find(c==c_v);
        
        % Check error.
        if r~=r_v(row_fit) || c~=c_v(col_fit)
            error('Check r(c) and r(c)_v(row(col)_fit).')
        end
        
        % Vectorized coordinate for fit.
        rr_v = kron(ones(length(r_v),1),r_v);
        cc_v = kron(c_v,ones(length(c_v),1));
        %ind_v = sub2ind(size(Smap_m),rr_v,cc_v);                
        
        % Get fit region data.
        s_m = Smap_m(r_v,c_v);
        h_m = I_body_m(r_v,c_v);
        bw_m = BWfit(r_v,c_v);
        mod_m = abs(h_m./s_m);
        mod_m(isinf(mod_m)) = 0;
        
        
        %-----------------------------------------
        % 2D fit.
        omega = n;
        w_m = exp(-( (r_v*ones(1,length(r_v)) - r).^2 + (ones(length(c_v),1)*c_v' - c).^2 )/omega^2) .* ...
            mod_m .* bw_m;
        
        breal_v = [];
        for l = 0:maxOrder
            for m = 0:maxOrder
                breal_v = [breal_v; sum(w_m(:).*real(s_m(:)).*(rr_v-r).^l .* (cc_v-c).^m)];                
            end
        end
        
        bimag_v = [];
        for l = 0:maxOrder
            for m = 0:maxOrder
                bimag_v = [bimag_v; sum(w_m(:).*imag(s_m(:)).*(rr_v-r).^l .* (cc_v-c).^m)];                
            end
        end
        
        A_m = zeros((maxOrder+1)^2,(maxOrder+1)^2);
        for l = 0:maxOrder
            for m = 0:maxOrder
                indr = (maxOrder+1)*l+m+1;
                for l1 = 0:maxOrder
                    for m1 = 0:maxOrder
                        indc = (maxOrder+1)*l1+m1+1;                        
                        a_v = [];
                        A_m(indr,indc) = [a_v, sum(w_m(:).*(rr_v-r).^(l+l1) .* (cc_v-c).^(m+m1))];                        
                    end
                end
            end
        end        
        
        % There is no difference between using QR and pinv.        
        %%% QR.
%         [Q,R] = qr(A_m);
%         preal = R \ (Q' * breal_v);
%         pimag = R \ (Q' * bimag_v);        
        %%% Peudoinverse.
        preal = pinv(A_m)*breal_v;
        pimag = pinv(A_m)*bimag_v;        
        
        % Get fitted value for voxel [r,c].
        val_v = [];
        for l=0:maxOrder
            for m=0:maxOrder
                val_v = [val_v, (rr_v-r).^l .* (cc_v-c).^m];
            end
        end
        z_m = reshape(val_v*preal + 1i*val_v*pimag,size(s_m));        
        z = z_m(row_fit,col_fit);
        %-----------------------------------------
        
                
        % Output for each coil.
        Smap_fit_m(r,c) = z;
        
    end % for ind_vox
    close(h)
    
    % Test show.
    %showimage(jet,abs(Smap_fit_m.*imerode(BWfit,strel('disk',3,0))),abs(Smap_fit_m.*BW))
    %showimage(jet,angle(Smap_fit_m.*imerode(BWfit,strel('disk',3,0))),angle(Smap_fit_m.*BW))
    
    % Output.
    Smap_fit_3d(:,:,coil) = Smap_fit_m;
    
    % Report.
    t = toc;
    fprintf('        finished coil[%d] in [%f]secs.\n',coil,t)    
  
end % for coil

% Clear.
clear  s_m  h_m  mod_m  z_m  Smap_fit_m  Smap_m  Smap_thresh_3d
clear  BW0  BW  BWrim  bw_m  w_m  Q  R  A*  b*  r*  c*  ind*



%% Generate noise correlation
% This is by (de Zwart et al., MRM 2002;48:1011 App.D)

% psi_m = zeros(Nc,'single');
% for indc1 = 1:Nc
%     for indc2 = 1:Nc
%         p1 = I_ref_3d(:,:,indc1);
%         p2 = I_ref_3d(:,:,indc2);        
%         p1 = p1(:);
%         p2 = p2(:);
%         
%         numerator = p1'*p2;
%         denominator = sqrt( (p1'*p1) * (p2'*p2) );
%         
%         psi_m(indc1,indc2) = numerator / denominator;
%     end
% end
% 
% % Clear.
% clear  p1  p2  numerator  denominator



%% Output
S = Smap_fit_3d;
% PSI = psi_m;
BW = BWfit;

clear  Smap_fit_3d  I_ref_3d  I_body_m  BWfit
%pack

fprintf('        S generated.\n')


%% END






