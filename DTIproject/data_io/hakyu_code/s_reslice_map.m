
%[s_reslice_map] reslices SMAP and FMAP to match slice positions of DWI
%scan.
%
%
% Last modified
% 2010.08.09
%   This incorporates,
%       getting corresponding slice positions of SMAP and FMAP to DWI
%       reslicing of SMAP and FMAP to DWI
% 2010.08.18.
%   Pack memory.
% 2010.09.22.
%   Add ~isempty(FMAPparams.filename).
% 2011.08.05.
%   Don't use slice left and right position for resampling. Just use the
%   central one.
% 2012.02.01.
%   Take care of scan mode = 3D case.
% 2012.05.23.
%   Take care of volume coil case (REFparams.filename is empty).
% 2012.06.01.
%   Take care of nan case for sl_smap_target_v and sl_fmap_target_v.
% 2012.06.16.
%   Take care of number of slices (kz) for 3D case.
%
% Ha-Kyu



%% Resample slices of SMAP and FMAP for NAV and IMG recon

fprintf('Resample slices of SMAP and FMAP for NAV and IMG recon\n')



%% Check if SMAP is already resliced for IMG and NAV

cd(saveDataDir_s)
if ~(exist('I_smap_img_4d.mat','file')==2) && ...
        ~(exist('I_smap_nav_4d.mat','file')==2) && ...
        ~isempty(REFparams.filename)
    load  I_smap_4d
    load  I_mask_S_3d
    [sny,snx,snz,snc] = size(I_smap_4d);
    
    
    % Get corresponding slice positions for SMAP
    fprintf('    Get corresponding SMAP slice positions\n')
    
    % Set corresponding DWI slices.
    if strcmpi(DWIparams.scan_mode,'3D')
        nSLICE = RECONparams.ENCima.oversample_resolutions(3) * DWIparams.chunks;
    else
        nSLICE = DWIparams.nSLICE;
    end
    sl_dwi_left = floor(nSLICE/2);
    sl_dwi_right = ceil(nSLICE/2)-1;
    sl_dwi_center_mm_v = linspace(-sl_dwi_left*(DWIparams.slice_thickness+DWIparams.slice_gap), ...
        sl_dwi_right*(DWIparams.slice_thickness+DWIparams.slice_gap),nSLICE);
    
    % SMAP slices.
    sl_smap_left = floor(snz/2); % including zero position
    sl_smap_right = ceil(snz/2)-1;
    sl_smap_center_mm_v = linspace(-sl_smap_left*REFparams.slice_thickness,...
        sl_smap_right*REFparams.slice_thickness,snz);
    sl_smap_v = 1:snz;
    sl_smap_target_v = interp1(sl_smap_center_mm_v,sl_smap_v,sl_dwi_center_mm_v);
    if length(sl_smap_target_v)~=nSLICE
        error('s_reslice_map:main','sl_smap_target_v or sl_smap_ref_v must equal to nSLICE.')
    end
    
    % Report.
    fprintf('      sl_smap_target_v\n')
    disp(sl_smap_target_v)
    
    % If NaN occured (when SMAP coverage is less than IMG).
    if any(isnan(sl_smap_target_v))
        ind_v = find(isnan(sl_smap_target_v));
        %* find from start
        for ind = 1:length(ind_v)
            a = sl_smap_target_v(ind_v(ind)+1);
            if ~isnan(a)
                val1 = a;
                n1 = ind_v(ind);
                break
            end
        end
        %* find from end
        for ind =length(ind_v):-1:1
            a = sl_smap_target_v(ind_v(ind)-1);
            if ~isnan(a)
                val2 = a;
                n2 = ind_v(ind);
                break
            end
        end
        sl_smap_target_v(1:n1) = val1;
        sl_smap_target_v(end:-1:n2) = val2;
    end
        
    
    % Resample SMAP slices
    fprintf('    Resample SMAP slices\n')
    
    [x,y,z] = meshgrid(1:snx,(1:sny)',1:snz);
    I_smap_resample_center_4d = zeros(sny,snx,length(sl_smap_target_v),snc,'single');
    I_smap_resample_left_4d = I_smap_resample_center_4d * 0;
    I_smap_resample_right_4d = I_smap_resample_center_4d * 0;
    
    for ind_coil = 1:snc
        %* Resample at the slice center.
        [xi,yi,zi] = meshgrid(1:snx,(1:sny)',sl_smap_target_v);
        I_smap_resample_center_4d(:,:,:,ind_coil) = ...
            interp3(x,y,z, I_smap_4d(:,:,:,ind_coil),xi,yi,zi);
        if ind_coil==1
            I_mask_resample_center_3d = interp3(x,y,z,I_mask_S_3d,xi,yi,zi);
        end
        
        %* Resample at the slice left.
%         sl_dwi_left_mm_v = sl_dwi_center_mm_v - DWIparams.slice_thickness/2;
%         sl_smap_target_left_v = interp1(sl_smap_center_mm_v,sl_smap_v,sl_dwi_left_mm_v);
%         [xi,yi,zi] = meshgrid(1:snx,(1:sny)',sl_smap_target_left_v);
%         I_smap_resample_left_4d(:,:,:,ind_coil) = ...
%             interp3(x,y,z, I_smap_4d(:,:,:,ind_coil),xi,yi,zi);
%         if ind_coil==1
%             I_mask_resample_left_3d = interp3(x,y,z,I_mask_S_3d,xi,yi,zi);
%         end
        
        %* Resample at the slice right.        
%         sl_dwi_right_mm_v = sl_dwi_center_mm_v + DWIparams.slice_thickness/2;
%         sl_smap_target_right_v = interp1(sl_smap_center_mm_v,sl_smap_v,sl_dwi_right_mm_v);
%         [xi,yi,zi] = meshgrid(1:snx,(1:sny)',sl_smap_target_right_v);
%         I_smap_resample_right_4d(:,:,:,ind_coil) = ...
%             interp3(x,y,z, I_smap_4d(:,:,:,ind_coil),xi,yi,zi);
%         if ind_coil==1
%             I_mask_resample_right_3d = interp3(x,y,z,I_mask_S_3d,xi,yi,zi);
%         end
    end
    
    %* Average all three resamples.
%     I_smap_resample_4d = (I_smap_resample_left_4d + I_smap_resample_center_4d + ...
%         I_smap_resample_right_4d) / 3;
%     I_mask_resample_3d = (I_mask_resample_left_3d + I_mask_resample_center_3d + ...
%         I_mask_resample_right_3d) / 3;
%     clear  I_smap_resample_left_4d  I_smap_resample_center_4d  I_smap_resample_right_4d
%     clear  I_mask_resample_left_3d  I_mask_resample_center_3d  I_mask_resample_right_3d
    
    %* Don't average three samples. Just use central one.
    I_smap_resample_4d = I_smap_resample_center_4d;
    I_mask_resample_3d = I_mask_resample_center_3d;  
    clear  I_smap_resample_center_4d  I_mask_resample_center_3d
    
    
    % Remove holes in mask.
    for ind_slice = 1:size(I_mask_resample_3d,3)
        I_mask_resample_3d(:,:,ind_slice) = imfill(I_mask_resample_3d(:,:,ind_slice),'holes');
    end
    I_mask_resample_3d = (I_mask_resample_3d > 0);
    
    % Reserve resampled SMAP and MASK.
    I_smap_nav_4d = I_smap_resample_4d;
    I_smap_img_4d = I_smap_resample_4d;
    I_mask_S_nav_3d = I_mask_resample_3d;
    I_mask_S_img_3d = I_mask_resample_3d;    
    clear  I_smap_resample_4d  I_mask_resample_3d  x  y  z  xi  yi  zi
    
    % Get new size.
    [sny,snx,snz,snc] = size(I_smap_img_4d);
    
    % Save S.
    cd(saveDataDir_s)
    eval(sprintf('save  I_smap_nav_4d    I_smap_nav_4d'))
    eval(sprintf('save  I_mask_S_nav_3d  I_mask_S_nav_3d'))
    eval(sprintf('save  I_smap_img_4d    I_smap_img_4d'))
    eval(sprintf('save  I_mask_S_img_3d  I_mask_S_img_3d'))
    
    % Report.
    fprintf('\n')
    fprintf('      SMAP to be used are resampled and saved\n')
    fprintf('        I_smap_nav_4d,     [sny,snx,snz,snc]  = [%d,%d,%d,%d]\n',size(I_smap_nav_4d))
    fprintf('        I_mask_S_nav_3d,   [sny,snx,snz]      = [%d,%d,%d]\n',size(I_mask_S_nav_3d))
    fprintf('        I_smap_img_4d,     [sny,snx,snz,snc]  = [%d,%d,%d,%d]\n',size(I_smap_img_4d))
    fprintf('        I_mask_S_img_3d,   [sny,snx,snz]      = [%d,%d,%d]\n',size(I_mask_S_img_3d))
    
    % Clear data after processing.
    clear  I_smap_*  I_mask_*
end
fprintf('\n')



%% Check if FMAP is already resliced for FMAP correction

cd(saveDataDir_s)
if ~(exist('I_fmap_4d.mat','file')==2) && ~isempty(FMAPparams.filename)
    load  I_fmap    
    I_fmap_4d = single(I_fmap.fmap);
    I_mask_F_3d = I_fmap.mask;
    [fny,fnx,fnz,fnt] = size(I_fmap_4d);
    
    
    % Get corresponding slice positions for FMAP
    fprintf('    Get corresponding FMAP slice positions\n')
    
    % Set corresponding slices.
    % Set corresponding DWI slices.
    if strcmpi(DWIparams.scan_mode,'3D')
        nSLICE = RECONparams.ENCima.oversample_resolutions(3) * DWIparams.chunks;
    else
        nSLICE = DWIparams.nSLICE;
    end
    sl_dwi_left = floor(nSLICE/2);
    sl_dwi_right = ceil(nSLICE/2)-1;
    sl_dwi_center_mm_v = linspace(-sl_dwi_left*(DWIparams.slice_thickness+DWIparams.slice_gap), ...
        sl_dwi_right*(DWIparams.slice_thickness+DWIparams.slice_gap),nSLICE);
    
    % FMAP slices.
    sl_fmap_left = floor(fnz/2); % including zero position
    sl_fmap_right = ceil(fnz/2)-1;
    sl_fmap_center_mm_v = linspace(-sl_fmap_left*FMAPparams.slice_thickness,...
        sl_fmap_right*FMAPparams.slice_thickness,fnz);
    sl_fmap_v = 1:fnz;
    sl_fmap_target_v = interp1(sl_fmap_center_mm_v,sl_fmap_v,sl_dwi_center_mm_v);
    if length(sl_fmap_target_v)~=nSLICE
        error('s_reslice_map:main','sl_fmap_target_v or sl_fmap_ref_v must equal to nSLICE.')
    end
    
    % Report.
    fprintf('      sl_fmap_target_v\n')
    disp(sl_fmap_target_v)
    
    % If NaN occured (when FMAP coverage is less than IMG).
    if any(isnan(sl_fmap_target_v))
        ind_v = find(isnan(sl_fmap_target_v));
        %* find from start
        for ind = 1:length(ind_v)
            a = sl_fmap_target_v(ind_v(ind)+1);
            if ~isnan(a)
                val1 = a;
                n1 = ind_v(ind);
                break
            end
        end
        %* find from end
        for ind =length(ind_v):-1:1
            a = sl_fmap_target_v(ind_v(ind)-1);
            if ~isnan(a)
                val2 = a;
                n2 = ind_v(ind);
                break
            end
        end
        sl_fmap_target_v(1:n1) = val1;
        sl_fmap_target_v(end:-1:n2) = val2;
    end
    
    
    % Resample FMAP slices.
    fprintf('    Resample FMAP slices\n')
    
    [x,y,z] = meshgrid(1:fnx,(1:fny)',1:fnz);
    I_fmap_resample_center_4d = zeros(fny,fnx,length(sl_fmap_target_v),fnt,'single');
    I_fmap_resample_left_4d = I_fmap_resample_center_4d * 0;
    I_fmap_resample_right_4d = I_fmap_resample_center_4d * 0;
    
    for ind_type = 1:fnt
        %* Resample at the slice center.
        [xi,yi,zi] = meshgrid(1:fnx,(1:fny)',sl_fmap_target_v);
        I_fmap_resample_center_4d(:,:,:,ind_type) = ...
            interp3(x,y,z, I_fmap_4d(:,:,:,ind_type),xi,yi,zi);
        I_mask_resample_center_3d = interp3(x,y,z,I_mask_F_3d,xi,yi,zi);
        
        %* Resample at the slice left.
%         sl_dwi_left_mm_v = sl_dwi_center_mm_v - DWIparams.slice_thickness/2;
%         sl_fmap_target_left_v = interp1(sl_fmap_center_mm_v,sl_fmap_v,sl_dwi_left_mm_v);
%         [xi,yi,zi] = meshgrid(1:fnx,(1:fny)',sl_fmap_target_left_v);
%         I_fmap_resample_left_4d(:,:,:,ind_type) = ...
%             interp3(x,y,z, I_fmap_4d(:,:,:,ind_type),xi,yi,zi);
%         I_mask_resample_left_3d = interp3(x,y,z,I_mask_F_3d,xi,yi,zi);
        
        %* Resample at the slice right.        
%         sl_dwi_right_mm_v = sl_dwi_center_mm_v + DWIparams.slice_thickness/2;
%         sl_fmap_target_right_v = interp1(sl_fmap_center_mm_v,sl_fmap_v,sl_dwi_right_mm_v);
%         [xi,yi,zi] = meshgrid(1:fnx,(1:fny)',sl_fmap_target_right_v);
%         I_fmap_resample_right_4d(:,:,:,ind_type) = ...
%             interp3(x,y,z, I_fmap_4d(:,:,:,ind_type),xi,yi,zi);
%         I_mask_resample_right_3d = interp3(x,y,z,I_mask_F_3d,xi,yi,zi);
    end
    
    %* Average all three resamples.
%     I_fmap_resample_4d = (I_fmap_resample_left_4d + I_fmap_resample_center_4d + ...
%         I_fmap_resample_right_4d) / 3;
%     I_mask_resample_3d = (I_mask_resample_left_3d + I_mask_resample_center_3d + ...
%         I_mask_resample_right_3d) / 3;
%     clear  I_fmap_resample_left_4d  I_fmap_resample_center_4d  I_fmap_resample_right_4d
%     clear  I_mask_resample_left_3d  I_mask_resample_center_3d  I_mask_resample_right_3d

    %* Don't average three samples. Just use central one.
    I_fmap_resample_4d = I_fmap_resample_center_4d;
    I_mask_resample_3d = I_mask_resample_center_3d;
    clear  I_fmap_resample_center_4d  I_mask_resample_center_3d
    
    
    % Remove holes in mask.
    for ind_slice = 1:size(I_mask_resample_3d,3)
        I_mask_resample_3d(:,:,ind_slice) = imfill(I_mask_resample_3d(:,:,ind_slice),'holes');
    end
    I_mask_resample_3d = (I_mask_resample_3d > 0);
    
    % Reserve resampled FMAP and MASK.
    I_fmap_4d = I_fmap_resample_4d;
    I_mask_F_3d = I_mask_resample_3d;    
    clear  I_fmap_resample_4d  I_mask_resample_3d  x  y  z  xi  yi  zi
    
    % Get new size.
    [fny,fnx,fnz,fnt] = size(I_fmap_4d);
    
    % Save S.
    cd(saveDataDir_s)
    eval(sprintf('save  I_fmap_4d    I_fmap_4d'))
    eval(sprintf('save  I_mask_F_3d    I_mask_F_3d'))
    
    % Report.
    fprintf('\n')
    fprintf('      FMAP to be used are resampled and saved\n')
    fprintf('        I_fmap_4d,       [fny,fnx,fnz,fnt]  = [%d,%d,%d,%d]\n',size(I_fmap_4d))
    fprintf('        I_mask_F_3d,     [fny,fnx,fnz]      = [%d,%d,%d]\n',size(I_mask_F_3d))
    
    % Clear data after processing.
    clear  I_fmap*  I_mask_*    
end



%% Clear data
clear  sl_smap_left  sl_smap_right  sl_fmap_left  sl_fmap_right
clear  sl_dw_left  sl_dwi_right  sl_dwi_v

% Pack memory.
pack

fprintf('\n\n')










