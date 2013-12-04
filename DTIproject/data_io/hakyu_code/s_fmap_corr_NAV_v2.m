
%[s_fmap_corr_NAV] does fieldmap correction of NAV to match IMG.
%
%
% Last modified
% 2010.02.11.
% 2010.03.23. Use resampled I_fmap_4d.
% 2010.07.28. Put (BW_img-BW_nav) / (BW_img*BW_nav) for doing
%             forward-backward distortion correction by one-time.
% 2010.08.02. Put (BW_img*BW_nav)/(BW_img-BW_nav) because this will be
%             divided in [fmap_corr_v4.m].
% 2010.08.09
%   This is generated from [NAV_fmapCorr.m].
%   The fmap correction will be done on windowed NAV recon image only,
%   because always the windowed version will be used for the final IMG
%   reconstruction with phase correction. Also fmap correction doesn't
%   require high-resolution (no windowed) NAV recon image data. Then only
%   windowed NAV recon images will be used here.
%   Use [f_fmap_corr.m] which is substitue for [fmap_corr_v4.m].
% 2010.08.18.
%   Pack memory.
% 2010.08.22.
%   Try different NAV, FMAP data with bandwidths. See note page 128.
%   Add flag_method for choosing between 'theoretical' and 'empirical'
%   bandwidth to be used for the correction.
% 2010.08.24.
%   Do FMAP correction to both of the non-windowed and windowed NAV recon
%   data. The reason for doing the correction for non-windowed NAV recon
%   data is so that the corrected non-windowed data can be used for the
%   recon with NAV image data.
% 2010.09.21.
%   Consider multiple non-zero b-values.
% 2010.10.27.
%   Change flag_method from 'empirical' to 'theoretical' as of 2010.10.26
%   for r2359 (3T) and r2336 (7T). This has not been tested yet but has
%   correct navigator acquisition bandwidth.
% 2010.12.07.
%   Modify 'blipPolarity' and 'slice_orientation'.
% 2011.01.04.
%   Add slice-by-slice co-registration between fmap_magn (fmap_freq) and
%   IMG (b=0) images.
% 2011.02.02.
%   Co-registration using my function, not vuTool. This is for
%   co-registeration between FMAP and IMG for FMAP correction of NAV.
% 2011.03.28.
%   A new method is added. But use old method with modified
%   FMAPparams.vox_kern in [s_get_FMAPparams.m].
% 2011.06.23.
%   Modified for transverse-RL case and others using img_reori_3d and
%   mask_reori.
% 2012.06.16.
%   Take care of scan_mode = 3D case.
% 2012.07.11.
%   NAV image is also coregistered to IMG. This is when multichunk (> 3) 3D
%   DTI acquisition is acquired. NAV is shifted alogn kx compared with IMG.
%   This must be corrected using some reference scan (FRC correction?) but
%   the data is not available in .LIST/.DATA data. So, as a temporary
%   solution, use coregistration.
%
%   Using this, recon_nav_ori##__R3S4 (e.g.) and recon_nav_win_ori##__R3S4
%   would be in native NAV space, but recon_nav_fmap_ori##__R3S4 and
%   recon_nav_win_fmap_ori##__R3S4 would be in coregistration with IMG,
%   only when scan_mode==3D.
% 2012.07.27.
%   Version 2 is generated.
%
%
% Ha-Kyu



%% Fieldmap correction of NAV for all nDW_GRAD and nSLICE
warning off

% Load fmap.
cd(saveDataDir_s)
load  I_fmap_4d
load  I_mask_F_3d
fprintf('\n')
fprintf('Load I_fmap_4d, I_mask_F_3d.\n')


% Load IMG (b=0) image to reorient it before coregistration between IMG and
% FMAP for NAV FMAP correction.
cd(saveReconDir_s)
load(sprintf('recon_img_ori00__R%dS%d',R,Ns));

for ind_chunk = 1:DWIparams.chunks
    
    if strcmpi(DWIparams.scan_mode,'3D')
        img_3d = img_recon_5d(:,:,:,ind_chunk,1); % [y,x,slice,chunk,avg]
    else
        img_3d = img_recon_5d(:,:,:,1); % [y,x,slice,avg]
    end
    
    
    % Get masked magnitude image.
    diskrad = 9; % disk radius
    mask = abs(img_3d)>0;
    [ny,nx,nz] = size(mask);
    for ind=1:nz
        mask(:,:,ind) = imerode(mask(:,:,ind),strel('disk',diskrad));
    end
    fprintf('    IMG(b=0) for chunk(%d) is loaded for co-registration between FMAP and IMG\n',...
        ind_chunk)
    cd(saveDataDir_s)
    
    
    % Revert object position of IMG done in [s_recon_IMG.m]. This is because
    % FMAP is in the same position as DWI data is reconstructed.
    if strcmpi(DWIparams.slice_orientation,'transverse')
        if strcmpi(DWIparams.fold_over_dir,'AP')
            switch DWIparams.fat_shift_dir
                case 'P'
                    % do nothing
                    img_reori_3d = img_3d;
                    mask_reori = mask;
                case 'A'
                    img_reori_3d = zeros(ny,nx,nz); % re-oriented image
                    mask_reori = zeros(ny,nx,nz);
                    for ind=1:nz
                        img_reori_3d(:,:,ind) = fliplr(flipud(img_3d(:,:,ind)));
                        mask_reori(:,:,ind) = fliplr(flipud(mask(:,:,ind)));
                    end
                otherwise
                    error('s_fmap_corr_NAV:main','Undefined fat_shift_dir')
            end
        elseif strcmpi(DWIparams.fold_over_dir,'RL')
            switch DWIparams.fat_shift_dir
                case 'L'
                    img_reori_3d = zeros(nx,ny,nz); % re-oriented image
                    mask_reori = zeros(nx,ny,nz);
                    for ind=1:nz
                        img_reori_3d(:,:,ind) = permute(flipud(img_3d(:,:,ind)),[2,1]);
                        mask_reori(:,:,ind) = permute(flipud(mask(:,:,ind)),[2,1]);
                    end
                case 'R'
                    img_reori_3d = zeros(nx,ny,nz); % re-oriented image
                    mask_reori = zeros(nx,ny,nz);
                    for ind=1:nz
                        img_reori_3d(:,:,ind) = permute(fliplr(img_3d(:,:,ind)),[2,1]);
                        mask_reori(:,:,ind) = permute(fliplr(mask(:,:,ind)),[2,1]);
                    end
                otherwise
                    error('s_fmap_cor_NAV:main','Undefined fat_shift_dir')
            end
        else
            error('s_fmap_corr_NAV:main','Undefined fold_over_dir')
        end
        
    elseif strcmpi(DWIparams.slice_orientation,'sagittal')
        
        if strcmpi(DWIparams.fold_over_dir,'AP')
            switch DWIparams.fat_shift_dir
                case 'P'
                    img_reori_3d = zeros(ny,nx,nz); % re-oriented image
                    mask_reori = zeros(ny,nx,nz);
                    for ind=1:nz
                        img_reori_3d(:,:,ind) = fliplr(img_3d(:,:,ind));
                        mask_reori(:,:,ind) = fliplr(mask(:,:,ind));
                    end
                case 'A'
                    img_reori_3d = zeros(ny,nx,nz); % re-oriented image
                    mask_reori = zeros(ny,nx,nz);
                    for ind=1:nz
                        img_reori_3d(:,:,ind) = fliplr(flipud(img_3d(:,:,ind)));
                        mask_reori(:,:,ind) = fliplr(flipud(mask(:,:,ind)));
                    end
                otherwise
                    error('s_fmap_corr_NAV:main','Undefined fat_shift_dir')
            end
        elseif strcmpi(DWIparams.fold_over_dir,'FH')
            warning('s_fmap_corr_NAV:main','SAGITTAL FH is not setup yet')
        else
            error('s_fmap_corr_NAV:main','Undefined fold_over_dir')
        end
        
    elseif strcmpi(DWIparams.slice_orientation,'coronal')
        
        if strcmpi(DWIparams.fold_over_dir,'RL')
            switch DWIparams.fat_shift_dir
                case 'L'
                    % do nothing
                    img_reori_3d = img_3d;
                    mask_reori = mask;
                case 'R'
                    img_reori_3d = img_3d;
                    mask_reori = mask;
                    warning('s_fmap_corr_NAV:main','CORONAL RL-R is not setup yet')
                otherwise
                    error('s_fmap_corr_NAV:main','Undefined fat_shift_dir')
            end
        elseif strcmpi(DWIparams.fold_over_dir,'FH')
            img_reori_3d = img_3d;
            mask_reori = mask;
            warning('s_fmap_corr_NAV:main','CORONAL FH is not setup yet')
        else
            error('s_fmap_corr_NAV:main','Undefined fold_over_dir')
        end
        
    else
        error('s_fmap_corr_NAV:main','Undefined slice_orientation')
    end
    
    img_3d = img_reori_3d;
    mask = mask_reori;
    clear  img_reori_3d  mask_reori
    
    % Make img_3d for each chunk.
    if strcmpi(DWIparams.scan_mode,'3D')
        eval(sprintf('img_chunk%.2d_3d = img_3d;',ind_chunk))
        eval(sprintf('mask_chunk%.2d_3d = mask;',ind_chunk))
        clear img_3d  mask
    end
    
end % for ind_chunk
clear  img_recon_5d


% FMAP correction of non-windowed and windowed NAV recon.
for dw_ori = 0:DWIparams.nDW_GRAD
    
    % Load NAV recon data without or with k-space window.
    if dw_ori==0
        win_v = 0:1;
        %win_v = 1;
    else
        win_v = 1;
    end
    for ind_nav_win = win_v
        if ind_nav_win==0
            % Load non-windowed NAV recon data.
            cd(saveReconDir_s)
            fname_s = sprintf('recon_nav_ori%.2d__R%dS%d',dw_ori,R,Ns);
            eval(sprintf('load  %s',fname_s))
            nav_recon_data = nav_recon_6d;
            clear  nav_recon_6d  nav_gfactor_6d
            fprintf('    [%s] loaded\n',fname_s)
        else
            % Load windowd NAV recon data.
            cd(saveReconDir_s)
            fname_s = sprintf('recon_nav_win_ori%.2d__R%dS%d',dw_ori,R,Ns);
            
            
            %----- TEMP -----
            flag = 1;
            while flag
                if ~exist(sprintf('%s.mat',fname_s),'file')
                    fprintf('    wait for 10 minutues\n')
                    pause(60*10)
                    flag = 1;
                else
                    s=dir(sprintf('%s.mat',fname_s));
                    if s.bytes > 300*10^6
                        fprintf('    do fmap corr\n')
                        flag = 0;
                    else
                        flag=1;
                    end
                end
            end
            %----- TEMP -----
            
            
            eval(sprintf('load  %s',fname_s))
            nav_recon_data = nav_recon_win_6d;
            clear  nav_recon_win_6d  nav_gfactor_win_6d
            fprintf('    [%s] loaded\n',fname_s)
        end
        if strcmpi(DWIparams.scan_mode,'3D')
            [nny,nnx,nnz,nnchunk,nnavg,nndw,nns] = size(nav_recon_data);
        else
            [nny,nnx,nnz,nnavg,nndw,nns] = size(nav_recon_data);
            nnchunk = 1;
        end
        
        % Remove nan.
        nav_recon_data(isnan(nav_recon_data)) = 0;
        nav_recon_data(isinf(nav_recon_data)) = 0;
        mask_nav = f_gen_mask2(nav_recon_data,0,'erode',6,'disk');
        fprintf('    NAV mask generated\n')
        
        
        % Preliminary.
        switch DWIparams.fold_over_dir
            case 'AP'
                PEdir = 1; % row direction
                if strcmpi(DWIparams.fat_shift_dir,'P')
                    blipPolarity = -1;  % in Matlab coordinate, desirable direction of deformation
                else
                    blipPolarity = 1;
                end
            case 'RL'
                %PEdir = 2; % col direction
                PEdir = 1; % row direction: phase-encoding is along Matlab row dir also
                if strcmpi(DWIparams.fat_shift_dir,'L')
                    blipPolarity = -1;  % in Matlab coordinate, desirable direction of deformation
                else
                    blipPolarity = 1;
                end
            otherwise
                error('s_fmap_corr_NAV:main','DWIparams.fold_over_dir must be ''AP'' or ''RL''.')
        end
        
        switch DWIparams.slice_orientation
            case 'transverse'
                slice_orientation = 1; % row direction
            case 'coronal'
                slice_orientation = 2; % col direction
            case 'sagittal'
                slice_orientation = 3; % slice direction
            otherwise
                error('s_fmap_corr_NAV:main','Unknown slice_orientation.')
        end
        
        % Get parameters.
        BW_img = DWIparams.BW(1);
        BW_nav = DWIparams.BW(2);
        
        % Reserve output.
        nav_recon_data_fmap = zeros(size(nav_recon_data),'single');
        
        % Loop.
        for ind_avg = 1:nnavg
            for ind_dw = 1:nndw
                for ind_chunk = 1:nnchunk
                    if strcmpi(DWIparams.scan_mode,'3D')
                        eval(sprintf('img_3d = img_chunk%.2d_3d;',ind_chunk))
                        eval(sprintf('mask = mask_chunk%.2d_3d;',ind_chunk))
                    end
                    
                    for ind_slice = 1:nnz
                        
                        % Generate FMAP slice. This is 1:nkz*nchunk, but
                        % ind_slice is 1:nkz. Then FMAP slice index must be
                        % redefined here.
                        %ind_slice1 = ind_slice + nnz * (ind_chunk - 1);
                        
                        % For [s_reslice_smap_v2.m] with unreversed SMAP slice
                        % order prescribed in [s_get_SMAPparams.m] and
                        % [f_coilsensemap.m]. 2012.07.10.
                        ind_slice1 = (nnchunk-ind_chunk)*nnz + ind_slice;
                        
                        % Read FMAP slice.
                        fmap_magn = I_fmap_4d(:,:,ind_slice1,1);
                        fmap_freq = I_fmap_4d(:,:,ind_slice1,2);
                        fmap_mask = I_mask_F_3d(:,:,ind_slice1);
                        
                        
                        % Dilate mask.
                        dilateval = 3; % voxel
                        fmap_mask = imdilate(fmap_mask,strel('disk',dilateval));
                        
                        
                        %-------------------------------------------------
                        % Do co-registration between FMAP and IMG (b=0).
                        img = img_3d(:,:,ind_slice).*mask(:,:,ind_slice);
                        
                        %* Use following set up as default.
                        optim_s = 'powell';
                        method_s = 'rigid';
                        cost_s = 'nmi';
                        p0_v = [0 0,0]'; % 2-D rigid body
                        nBin = 256;
                        flag_display = 'off'; % 'on' or 'off'
                        
                        %* Co-registration.
                        
                        %* Co-registration only need to be done once for
                        %* dw_ori==0 data.
                        
                        if dw_ori==0 && ind_nav_win==0
                            tic
                            p1_v = f_powell_search('f_optim_affine_coreg',p0_v, ...
                                abs(img),fmap_magn,method_s,cost_s, ...
                                nBin,flag_display);                            
                            fprintf('    FMAP is co-registered to IMG using [%s]-[%s]-[%s]: ori[%d]-windowed[%d]-avg[%d]-dw[%d]-chunk[%d]-sl[%d] in [%.3f]sec\n',...
                                upper(optim_s),upper(method_s),upper(cost_s), ...
                                dw_ori,ind_nav_win,ind_avg,ind_dw,ind_chunk,ind_slice,toc)
                            eval(sprintf('save  p1_v__sl%.2d_chunk%.2d_dw%.2d_avg%.2d  p1_v',...
                                ind_slice,ind_chunk,ind_dw,ind_avg))
                            fprintf('      p1_v is saved for sl%.2d_chunk%.2d_dw%.2d_avg%.2d\n',...
                                ind_slice,ind_chunk,ind_dw,ind_avg)
                        else
                            eval(sprintf('load  p1_v__sl%.2d_chunk%.2d_dw%.2d_avg%.2d  p1_v',...
                                ind_slice,ind_chunk,ind_dw,ind_avg))
                            fprintf('      p1_v is loaded for sl%.2d,chunk%.2d,dw%.2d,avg%.2d,win%.2d,ori%.2d\n',...
                                ind_slice,ind_chunk,ind_dw,ind_avg,ind_nav_win,dw_ori)
                        end
                        
                        %* Generate co-registered image.
                        [x_m,y_m] = meshgrid(1:nnx,1:nny);
                        one_v = ones(1,nnx*nny);
                        C = [x_m(:)'; y_m(:)'; one_v];
                        
                        T = f_tform_affine(fmap_magn,p1_v,method_s);
                        Cp = T*C;
                        xp_m = reshape(Cp(1,:),nny,nnx);
                        yp_m = reshape(Cp(2,:),nny,nnx);
                        
                        fmap_magn_reg = interp2(x_m,y_m,fmap_magn,xp_m,yp_m);
                        fmap_magn_reg(isnan(fmap_magn_reg)) = 0;
                        fmap_magn_reg(isinf(fmap_magn_reg)) = 0;
                        fmap_freq_reg = interp2(x_m,y_m,fmap_freq,xp_m,yp_m);
                        fmap_freq_reg(isnan(fmap_freq_reg)) = 0;
                        fmap_freq_reg(isinf(fmap_freq_reg)) = 0;
                        
                        %* Use co-registered image for FMAP correction.
                        fmap_magn = fmap_magn_reg;
                        fmap_freq = fmap_freq_reg;
                        clear  fmap_magn_reg  fmap_freq_reg  p1_v
                        
                        %* Check.
                        %m1 = gen_contour_overlaid(abs(img),edge(fmap_magn,'canny'),0.3);
                        %m2 = gen_contour_overlaid(abs(img),edge(fmap_magn_reg,'canny'),0.3);
                        %showimage(jet,m1,m2,'before','after')
                        %-----------------END------------------------------
                        
                        
                        
                        %-------------------------------------------------
                        % Do co-registration between NAV and IMG (b=0).
                                                
                        %* Sometimes IMG and NAV are in different position
                        %* when 3D multichunk acquisition (chunk > 3) is
                        %* used. Then do coregistration between NAV and IMG
                        %* before FMAP correction of NAV.
                        
                        %--- DON'T DO THIS. THIS IS CAUSED BY WRONG NAV
                        %--- DATA SHIFT!
                        
                        if strcmpi(DWIparams.scan_mode,'3D') && ...
                                (DATAFLAGparams.coreg_nav_to_img==true)
                            
                            nav = squeeze(nav_recon_data(:,:,ind_slice, ...
                                ind_chunk,ind_avg,ind_dw,:) .* ...
                                mask_nav(:,:,ind_slice,ind_chunk,ind_avg,...
                                ind_dw,:));
                            
                            % Without mask, registration isn't good.
                            %nav = squeeze(nav_recon_data(:,:,ind_slice, ...
                            %    ind_chunk,ind_avg,ind_dw,:));
                            
                            %* Use following set up as default.
                            optim_s = 'powell';
                            method_s = 'rigid';
                            cost_s = 'nmi';
                            p0_v = [0 0,0]'; % 2-D rigid body
                            nBin = 128;
                            flag_display = 'off'; % 'on' or 'off'
                            
                            %* Co-registration.
                            for ind_shot = 1:size(nav,3)
                                tic
                                p1_v = f_powell_search('f_optim_affine_coreg',p0_v, ...
                                    abs(img),abs(nav(:,:,ind_shot)),method_s,cost_s, ...
                                    nBin,flag_display);
                                fprintf('      NAV is co-registered to IMG using [%s]-[%s]-[%s]: ori[%d]-windowed[%d]-avg[%d]-dw[%d]-chunk[%d]-sl[%d]-shot[%d] in [%.3f]sec\n',...
                                    upper(optim_s),upper(method_s),upper(cost_s), ...
                                    dw_ori,ind_nav_win,ind_avg,ind_dw, ...
                                    ind_chunk,ind_slice,ind_shot,toc)
                                
                                %* Generate co-registered image.
                                [x_m,y_m] = meshgrid(1:nnx,1:nny);
                                one_v = ones(1,nnx*nny);
                                C = [x_m(:)'; y_m(:)'; one_v];
                                
                                T = f_tform_affine(nav(:,:,ind_shot),p1_v,method_s);
                                Cp = T*C;
                                xp_m = reshape(Cp(1,:),nny,nnx);
                                yp_m = reshape(Cp(2,:),nny,nnx);
                                
                                nav_reg = interp2(x_m,y_m,nav(:,:,ind_shot),xp_m,yp_m);
                                nav_reg(isnan(nav_reg)) = 0;
                                nav_reg(isinf(nav_reg)) = 0;
                                
                                nav(:,:,ind_shot) = nav_reg;
                            end
                            clear  nav_reg
                            
                            % * Keep registered data.
                            nav_recon_data(:,:,ind_slice,ind_chunk,ind_avg,ind_dw,:) = ...
                                nav;
                            
                            
                            %* Check.
                            %m1 = gen_contour_overlaid(abs(img),edge(abs(nav(:,:,ind_shot)),'canny'),0.2);
                            %m2 = gen_contour_overlaid(abs(img),edge(abs(nav_magn_reg),'canny'),0.2);
                            %showimage(jet,m1,m2,'before','after')
                        end
                        %-----------------END------------------------------
                        
                        
                        % Smooth FMAP.
                        if strcmpi(FMAPparams.fitregion,'local')
                            %*** Local polynomial fit.
                            [F,BW] = polyfit_fieldmap(fmap_magn, fmap_freq, fmap_mask, ...
                                FMAPparams.fitorder, FMAPparams.fitsize(1), ...
                                FMAPparams.fitsize(2), FMAPparams.vox_dilate);
                            fmap_freq_fit = F;
                            fmap_mask = BW;
                            fmap_freq_fit = fmap_freq_fit.*fmap_mask;
                            clear  F  BW
                        elseif strcmpi(FMAPparams.fitregion,'global')
                            %*** Global polynomial fit.
                            kern_vox = FMAPparams.vox_kern;
                            kern = smoothkern([kern_vox,kern_vox],kern_vox/2);
                            fmap_freq_smooth = conv2(fmap_freq,kern,'same');
                            
                            se = strel('disk',FMAPparams.vox_dilate,0);    % vox_dilate for fieldmap mask
                            fmap_mask_dil = imdilate(fmap_mask,se);
                            fmap_freq_mask = fmap_freq_smooth.*fmap_mask_dil;
                            
                            maxOrder = FMAPparams.fitorder;
                            [r_v,c_v] = find(fmap_mask_dil);
                            ind_v = find(fmap_mask_dil(:));
                            pos_m = [];
                            for l=0:maxOrder
                                for n=0:maxOrder
                                    pos_m = [pos_m, (r_v).^l .* (c_v).^n];
                                end
                            end
                            p = pinv(pos_m)*fmap_freq_mask(ind_v);
                            fmap_freq_fit = zeros(size(fmap_freq_mask));
                            fmap_freq_fit(ind_v) = pos_m*p;
                            clear  fmap_mask  fmap_freq_mask  fmap_freq_smooth  fmap_mask_dil
                        elseif strcmpi(FMAPparams.fitregion,'none')
                            %*** No fit, but smooth.
                            
                            %------------ old method --------------
                            kern_vox = FMAPparams.vox_kern;
                            kern = f_gauss_kernel([kern_vox,kern_vox],kern_vox/2);
                            fmap_freq_smooth = conv2(fmap_freq,kern,'same');
                            %fmap_freq_fit = fmap_freq_smooth.*fmap_mask;
                            %fmap_freq_smooth = medfilt2(fmap_freq,[2,2]); % not better than using [f_gauss_kernel.m] above
                            fmap_freq_fit = fmap_freq_smooth;
                            
                            %------------ new method --------------
%                             v = fmap_freq(fmap_freq~=0);
%                             [n_v,x_v]=hist(v);
%                             vox_disp_max=x_v(n_v==max(n_v));
%                             kern_vox = ceil(abs(vox_disp_max/BW_nav));
%                             if kern_vox < 2 % min kernel size
%                                 kern_vox=2;
%                             elseif kern_vox > 4 % max kernel size
%                                 kern_vox=4;
%                             end
%                             kern = f_gauss_kernel([kern_vox,kern_vox],kern_vox/2);
%                             fmap_freq_smooth = conv2(fmap_freq,kern,'same');
%                             fmap_freq_fit = fmap_freq_smooth.*fmap_mask;
                        else
                            error('s_fmap_corr_NAV:main','    Unknown FMAPparams.fitregion.')
                        end
                        
                        
                        % Generate 3D fieldmap data, F.
                        F = squeeze(cat(3,fmap_magn,fmap_freq_fit));
                        clear  fmap_magn  fmap_freq_fit
                        
                        
                        %----------------------------------------------------------
                        % Fieldmap correction: forward = NAV space to original space
                        %blipDir = blipPolarity;
                        %nav_recon_fmap_fw_3d = zeros(size(squeeze(nav_recon_data(:,:,ind_slice,ind_avg,:))),'single');
                        %for ind_shot = 1:Ns
                        %    corr_data = f_fmap_corr(squeeze(nav_recon_data(:,:,ind_slice,ind_avg,ind_shot)), ...
                        %        F, BW_nav, PEdir, blipDir, slice_orientation);
                        %    nav_recon_fmap_fw_3d(:,:,ind_shot) = corr_data.Icorr;
                        %    clear  corr_data
                        %end
                        
                        
                        %----------------------------------------------------------
                        % Fieldmap correction: backward = original space to IMG space
                        %blipDir = -blipPolarity;
                        %nav_recon_fmap_fwbw_3d = zeros(size(nav_recon_fmap_fw_3d),'single');
                        %for ind_shot = 1:Ns
                        %    corr_data = f_fmap_corr(nav_recon_fmap_fw_3d(:,:,ind_shot), ...
                        %        F, BW_img, PEdir,blipDir,slice_orientation);
                        %    nav_recon_fmap_fwbw_3d(:,:,ind_shot) = corr_data.Icorr;
                        %    clear  corr_data
                        %end
                        
                        
                        %-----------------------------------
                        % Fieldmap correction: forward-backward, this can replace above
                        % two processes.
                        blipDir = blipPolarity;
                        if strcmpi(DWIparams.scan_mode,'3D')
                            nav_recon_fmap_fwbw_3d = ...
                                zeros(size(squeeze( ...
                                nav_recon_data(:,:,ind_slice,ind_chunk,ind_avg,ind_dw,:))),'single');
                        else
                            nav_recon_fmap_fwbw_3d = ...
                                zeros(size(squeeze( ...
                                nav_recon_data(:,:,ind_slice,ind_avg,ind_dw,:))),'single');
                        end
                        
                        
                        % Change flag_method from 'empirical' to 'theoretical'
                        % as of 2010.10.26 for r2359 (3T) and r2336 (7T). This
                        % has not been tested yet but has correct navigator
                        % acquisition bandwidth.
                        %
                        % This has been tested.
                        
                        flag_method = 'empirical'; % ['theoretical','empirical']
                        %flag_method = 'theoretical';
                        for ind_shot = 1:Ns                            
                            if strcmpi(flag_method,'theoretical')
                                % Theoretically correct routine.
                                if strcmpi(DWIparams.scan_mode,'3D')
                                    corr_data = f_fmap_corr(squeeze( ...
                                        nav_recon_data(:,:,ind_slice,ind_chunk,ind_avg,ind_dw,ind_shot)), ...
                                        F, (BW_img*BW_nav)/(BW_img-BW_nav), PEdir, blipDir, ...
                                        slice_orientation);
                                else
                                    corr_data = f_fmap_corr(squeeze( ...
                                        nav_recon_data(:,:,ind_slice,ind_avg,ind_dw,ind_shot)), ...
                                        F, (BW_img*BW_nav)/(BW_img-BW_nav), PEdir, blipDir, ...
                                        slice_orientation);
                                end
                                
                            elseif strcmpi(flag_method,'empirical')
                                % Empirical routine (before 7T r2336).
                                % Use BW_nav/4 as bandwidth as empirical bandwidth
                                % difference.
                                %corr_data = f_fmap_corr(squeeze( ...
                                %    nav_recon_data(:,:,ind_slice,ind_avg,ind_dw,ind_shot)), ...
                                %    F, BW_nav/4, PEdir, blipDir, ...
                                %    slice_orientation);
                                
                                % Empirical routine (after 7T r2336 including it).
                                %corr_data = f_fmap_corr(squeeze( ...
                                %    nav_recon_data(:,:,ind_slice,ind_avg,ind_dw,ind_shot)), ...
                                %    F, (BW_img*BW_nav)/(BW_img-BW_nav) / (DWIparams.nSHOT/4), PEdir, blipDir, ...
                                %    slice_orientation);
                                
                                % Multiply interp_factor along PE direction so
                                % that a little bit less correction is made.
                                if strcmpi(DWIparams.scan_mode,'3D')
                                    corr_data = f_fmap_corr(squeeze( ...
                                        nav_recon_data(:,:,ind_slice,ind_chunk,ind_avg,ind_dw,ind_shot)), ...
                                        F, (BW_img*BW_nav)/(BW_img-BW_nav) * RECONparams.ENCima.interp_factors(2), ...
                                        PEdir, blipDir, ...
                                        slice_orientation);
                                else
                                    corr_data = f_fmap_corr(squeeze( ...
                                        nav_recon_data(:,:,ind_slice,ind_avg,ind_dw,ind_shot)), ...
                                        F, (BW_img*BW_nav)/(BW_img-BW_nav) * RECONparams.ENCima.interp_factors(2), ...
                                        PEdir, blipDir, ...
                                        slice_orientation);
                                end
                                
                            else
                                error('s_fmap_corr_NAV:main','Unknown flag_method. It must be ''theoretical'' or ''empirical''.')
                            end
                            
                            nav_recon_fmap_fwbw_3d(:,:,ind_shot) = corr_data.Icorr;
                            clear  corr_data
                        end
                        
                        
                        %-----------------------------------
                        % Test show NAV shot1 data and fmap GRE data.
                        %bwe = edge(F(:,:,1),'canny');
                        %figure
                        %showimage(jet,circshift(F(:,:,1),[yshift_vox,0]), ...
                        %    abs(nav_recon_3d(:,:,1))+bwe*10, ...
                        %    abs(nav_recon_fmap_fw_3d(:,:,1))+bwe*10, ...
                        %    abs(nav_recon_fmap_fwbw_3d(:,:,1))+bwe*10,...
                        %    'Fmap','NAV recon shot1','NAV recon fmap(forward)', ...
                        %    'NAV recon fmap(forward-backward)')
                        
                        
                        %-----------------------------------
                        % Output.
                        if strcmpi(DWIparams.scan_mode,'3D')
                            nav_recon_data_fmap(:,:,ind_slice,ind_chunk,ind_avg,ind_dw,:) = nav_recon_fmap_fwbw_3d;
                        else
                            nav_recon_data_fmap(:,:,ind_slice,ind_avg,ind_dw,:) = nav_recon_fmap_fwbw_3d;
                        end
                        
                    end % for ind_slice
                end % for ind_chunk
            end % for ind_dw
        end % for ind_avg        
        
        
        % Save FMAP corrected NAV recon data without or with k-space
        % window.
        if ind_nav_win==0
            % Save non-windowed NAV data.
            cd(saveReconDir_s)
            nav_recon_fmap_6d = nav_recon_data_fmap;
            fname_s = sprintf('recon_nav_fmap_ori%.2d__R%dS%d',dw_ori,R,Ns);
            eval(sprintf('save  %s  nav_recon_fmap_6d',fname_s))
            fprintf('    [%s] is saved at \n',fname_s)
            fprintf('    [%s]]\n',saveReconDir_s)
            clear  nav_recon_*
        else
            % Save windowed NAV data.
            cd(saveReconDir_s)
            nav_recon_win_fmap_6d = nav_recon_data_fmap;
            fname_s = sprintf('recon_nav_win_fmap_ori%.2d__R%dS%d',dw_ori,R,Ns);
            eval(sprintf('save  %s  nav_recon_win_fmap_6d',fname_s))
            fprintf('    [%s] is saved at \n',fname_s)
            fprintf('    [%s]]\n',saveReconDir_s)
            clear  nav_recon_*
        end
        
    end % for ind_nav_win
    
end % for dw_ori


% Clear data.
eval(sprintf('clear  img_chunk%.2d_3d  mask_chunk%.2d_3d', ...
    ind_chunk,ind_chunk))
clear  I_fmap*  I_mask*  FOV  sense_fac_p  ovs_fac_p  acq_vox_p
clear  Tro  PEdir  blipDir  slice_orientation  F  fmap*  kern*  mask  img_3d


% Pack memory.
pack

fprintf('\n\n')



%% END






