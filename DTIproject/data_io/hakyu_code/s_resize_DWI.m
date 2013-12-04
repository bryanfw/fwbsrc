
%[s_resize_DWI] resizes the final recon DWI images to have the recon voxel
%size.
%
%
% Last modified
% 2010.08.10
% 2010.08.18.
%   Pack memory.
% 2010.08.20.
%   Check if resized dw_data_*_5d is already exist. -> DON'T USE THE
%   CHECKING.
% 2010.09.01.
%   Selectively take DWIparams.FOV.[AP,RL] based on fold-over direction.
% 2010.09.11.
%   cut_vox_post_P_ORI+N_recon_P_ORI is modified to
%   cut_vox_pre_P_ORI+N_recon_P_ORI.
% 2010.09.14.
%   Check if resized data is already exist.
% 2010.09.29.
%   navg is changed to ndw in generating dw_data_5d.
% 2010.10.11.
%   Resize for no phase correction data.
% 2010.11.11.
%   Add coronal slice.
% 2010.12.08.
%   Modify resizing method 1 and 2 (flag_resMethod). Then flag_resMethod of
%   1 and 2 generates the same results.
% 2010.12.22.
%   Apply round() to N_recon_P_ORI and N_recon_M_ORI for flag_resMethod=1.
% 2011.01.19.
%   Re-orienta resized images based on <default display orientation> as
%   shown in E4p44.
% 2011.10.19.
%   Modify fname_s for loading reconstructed data. This is because
%   [s_recon_IMG.m] is no longer recon all non-corrected data.
%   Non-corrected data is reconstructed in [s_recon_DWI.m] as corrected
%   data does.
% 2012.02.13.
%   Add flag for resize to recon resolution or not.
% 2012.03.10.
%   'flag_average_idx' is added to selectively resize for corresponding
%   average index. If it is [], then use all averages.
%   This is applied only for phaseCorr==0,resize_DWI==1 or
%   phaseCorr==1,navFmapCorr==1,resize_DWI==1 cases.
% 2012.05.24.
%   Take care of Volume Coil case, REFparams.filename is empty.
%
% Ha-Kyu



%% Check if resized data is already exist
cd(saveReconDir_s)
flag_resizeDWI = [];

% Check flag_average_idx.
if isempty(flag_average_idx)
    tag_average_idx = [];
else
    tag_average_idx = sprintf('_avgIdx%d',flag_average_idx);
end

% Check if data is already exist.
if RECONFLAGparams.phaseCorr==0
    if RECONFLAGparams.resize_DWI==1
        if exist(sprintf('dw_data_no_corr_5d%s.mat',tag_average_idx),'file')==2
            fprintf('[dw_data_no_corr_5d%s.mat] exists. Then skip resizing.\n',tag_average_idx)
            fprintf('\n\n')
            return
        else
            fprintf('Resize to generate [dw_data_no_corr_5d%s.mat]\n',tag_average_idx)
        end
    else
        if exist('dw_data_no_corr_no_res_5d.mat','file')==2
            fprintf('[dw_data_no_corr_no_res_5d.mat] exists. Then skip generating.\n')
            fprintf('\n\n')
            return
        else
            fprintf('Do not resize, just generate [dw_data_no_corr_no_res_5d.mat]\n')
        end
    end
end
if RECONFLAGparams.phaseCorr==1
    if RECONFLAGparams.navFmapCorr==0 && RECONFLAGparams.navDeform==0
        if RECONFLAGparams.resize_DWI==1
            if exist('dw_data_5d.mat','file')==2
                fprintf('[dw_data_5d.mat] exists. Then skip resizing.\n')
                fprintf('\n\n')
                return
            else
                fprintf('Resize to generate [dw_data_5d.mat]\n')
            end
        else
            if exist('dw_data_no_res_5d.mat','file')==2
                fprintf('[dw_data_no_res_5d.mat] exists. Then skip generating.\n')
                fprintf('\n\n')
                return
            else
                fprintf('Do no resize, just generate [dw_data_no_res_5d.mat]\n')
            end
        end
    end
    if RECONFLAGparams.navFmapCorr==1
        if RECONFLAGparams.resize_DWI==1
            if exist(sprintf('dw_data_fmap_5d%s.mat',tag_average_idx),'file')==2
                fprintf('[dw_data_fmap_5d%s.mat] exists. Then skip resizing.\n',tag_average_idx)
                fprintf('\n\n')
                return
            else
                fprintf('Resize to generate [dw_data_fmap_5d%s.mat]\n',tag_average_idx)
            end
        else
            if exist('dw_data_fmap_no_res_5d.mat','file')==2
                fprintf('[dw_data_fmap_no_res_5d.mat] exists. Then skip generating.\n')
                fprintf('\n\n')
                return
            else
                fprintf('Do not resize, just generate [dw_data_fmap_no_res_5d.mat]\n')
            end
        end
    end
    if RECONFLAGparams.navDeform==1
        if RECONFLAGparams.resize_DWI==1
            if exist('dw_data_deform_5d.mat','file')==2
                fprintf('[dw_data_deform_5d.mat] exists. Then skip resizing.\n')
                fprintf('\n\n')
                return
            else
                fprintf('Resize to generate [dw_data_deform_5d.mat]\n')
            end
        else
            if exist('dw_data_deform_no_res_5d.mat','file')==2
                fprintf('[dw_data_deform_no_res_5d.mat] exists. Then skip generating.\n')
                fprintf('\n\n')
                return
            else
                fprintf('Do not resize, just generate [dw_data_deform_no_res_5d.mat]\n')
            end
        end
    end
end



%% Determine which NAV recon image to use

% Determine which NAV recon image to use.
if RECONFLAGparams.phaseCorr == 0
    fname_nav_head_s = sprintf('recon_nav_win');
    dataname_nav_s = sprintf('nav_recon_win_6d');
    fname_img_head_s = sprintf('recon_img');
else
    fname_nav_head_s = sprintf('recon_nav_win');
    dataname_nav_s = sprintf('nav_recon_win_6d');
    fname_img_head_s = sprintf('recon_img_corr');
    
    if RECONFLAGparams.navFmapCorr==1
        fname_nav_head_s = sprintf('recon_nav_win_fmap');
        dataname_nav_s = sprintf('nav_recon_win_fmap_6d');
        fname_img_head_s = sprintf('recon_img_corr_fmap');
    end
    
    if RECONFLAGparams.navDeform==1
        fname_nav_head_s = sprintf('recon_nav_deform_win');
        dataname_nav_s = sprintf('nav_recon_deform_win_6d');
        fname_img_head_s = sprintf('recon_img_corr_deform');
    end
    
end

fprintf('NAV data used\n')
fprintf('  fname_nav_head_s: %s\n',fname_nav_head_s)
fprintf('  dataname_nav_s: %s\n',dataname_nav_s)
fprintf('  fname_img_head_s: %s\n',fname_img_head_s)
fprintf('\n')



%% Resize DWI recon data to isotropic in-plane voxel size

% Other params.
R = DWIparams.SENSE_FACTOR;
Ns = DWIparams.nSHOT;


% Load mask data.
cd(saveDataDir_s)
%if ~isempty(REFparams.filename)
if exist('I_smap_img_4d.mat','file')    
    load  I_mask_S_img_3d
    [Ny,Nx,Nz] = size(I_mask_S_img_3d);
    clear  I_mask_S_img_3d
end

% Take care of Volume Coil case (REFparams.filename is empty).
cd(saveReconDir_s)
if isempty(REFparams.filename)    
    load(sprintf('recon_img_ori00__R%dS%d.mat',R,Ns))
    [Ny,Nx,Nz,Navg] = size(img_recon_5d);
    clear  img_recon_5d
    I_smap_img_4d = ones(Ny,Nx,Nz,Nc);
    I_mask_S_img_3d = ones(Ny,Nx,Nz);
end

% Reserve output.
if RECONFLAGparams.resize_DWI==1
    ny = RECONparams.a_N_final_interp_recon_p;
    nx = ny;
    %dw_data_5d = zeros(ny,nx,Nz,1,DWIparams.nDW_GRAD+1,'single'); % navg is always one
    if strcmpi(DWIparams.scan_mode,'3D')
        nchunk = DWIparams.chunks;
        %nz = Nz/nchunk;
        nz = DWIparams.kz;
        dw_data_5d = zeros(ny,nx,nz,nchunk,DWIparams.nROW-1,DWIparams.nDW_GRAD+1,'single');
    else
        dw_data_5d = zeros(ny,nx,Nz,DWIparams.nROW-1,DWIparams.nDW_GRAD+1,'single'); % navg is changed to ndw
        nchunk = 1;
    end
else
    if strcmpi(DWIparams.scan_mode,'3D')
        nchunk = DWIparams.chunks;
        %nz = Nz/nchunk;
        nz = DWIparams.kz;
        dw_data_5d = zeros(Ny,Nx,nz,nchunk,DWIparams.nROW-1,DWIparams.nDW_GRAD+1,'single');
    else        
        dw_data_5d = zeros(Ny,Nx,Nz,DWIparams.nROW-1,DWIparams.nDW_GRAD+1,'single');
        nchunk = 1;
    end
end


% Load, resize, zerofill recon data.
for dw_ori = 0:DWIparams.nDW_GRAD
    
    % Load DWI recon data.
    
    % Note:
    % Complex average of b=0 images over avg looks the best compared with
    % other averages such as averaging magn and phase separately or even
    % recontruction with phase correction over all averages. The complex
    % averaged b=0 data is the same as the data reconstructed without phase
    % correction.
    
    cd(saveReconDir_s)
    if RECONFLAGparams.phaseCorr == 0 % no correction
        %         fname_s = sprintf('%s_ori%.2d__R%dS%d', ...
        %             fname_img_head_s,dw_ori,R,Ns);
        
        % Non-corrected data is now reconstructed in [s_recon_DWI.m] not in
        % [s_recon_IMG.m] for DW data. Then here is the new fname_s.
        if isempty(flag_average_idx)
            fname_s = sprintf('%s_avg%.2d_ori%.2d__R%dS%d', ...
                fname_img_head_s,DWIparams.nNSA,dw_ori,R,Ns);
        else
            fname_s = sprintf('%s_avg%.2d_ori%.2d__R%dS%d%s', ...
                fname_img_head_s,1,dw_ori,R,Ns,tag_average_idx);
        end
    elseif RECONFLAGparams.phaseCorr == 1 % correction
        if isempty(flag_average_idx)
            fname_s = sprintf('%s_avg%.2d_ori%.2d__R%dS%d', ...
                fname_img_head_s,DWIparams.nNSA,dw_ori,R,Ns);
        else
            fname_s = sprintf('%s_avg%.2d_ori%.2d__R%dS%d%s', ...
                fname_img_head_s,1,dw_ori,R,Ns,tag_average_idx);
        end
    else
        error('s_resize_DWI:main','Unknown ''RECONFLAGparams.phaseCorr''.')
    end
    load(fname_s)
    fprintf('Loading [%s]\n',fname_s)
    
    
    % Generate mask.
    %if dw_ori==0
%     if RECONFLAGparams.phaseCorr==0
%         [n1,n2,n3,n4,n5] = size(img_recon_5d); % [ny,nx,nz,navg,ndw]
%         
%         % check if ndw is 1.
%         if n5 > 1
%             fprintf('\n\n\n')
%             error('n5 is more than 1 at RECONFLAGparams.phaseCorr==0')
%         elseif n5 == 1
%             img_recon_4d = squeeze(img_recon_5d);
%         end
%     end
        %mask = f_gen_mask2(img_recon_4d,0,'erode',6,'disk');
        mask = f_gen_mask2(img_recon_4d,0,'erode',12,'disk');
    %end
    
    % TEMP-----
    %mask = f_gen_mask2(img_recon_4d,0,'erode',6,'disk');
    % TEMP-----
    
    % Get recon image size.
    if strcmpi(DWIparams.scan_mode,'3D')
        [ny,nx,nz,nchunk,ndw] = size(img_recon_4d);
    else
        [ny,nx,nz,ndw] = size(img_recon_4d);
        nchunk = 1;
    end
    
    
    % Resize recon DWI image to have recon voxel size.
    flag_resMethod = 1; % set resize method
    
    for ind_dw = 1:ndw
        for ind_chunk = 1:nchunk            
        for ind_slice = 1:nz
            
            % Resize or not.
            if RECONFLAGparams.resize_DWI==0
                %fprintf('  Do not resize image to recon voxel size\n')
                if strcmpi(DWIparams.scan_mode,'3D')
                    img_recon_m = img_recon_4d(:,:,ind_slice,ind_chunk,ind_dw) .* mask(:,:,ind_slice,ind_chunk,ind_dw);
                    dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = img_recon_m;
                else
                    img_recon_m = img_recon_4d(:,:,ind_slice,ind_dw) .* mask(:,:,ind_slice);
                    dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = img_recon_m;
                end
            else                
                %fprintf('  Resize image to recon voxel size\n')
                if strcmpi(DWIparams.scan_mode,'3D')
                    img_recon_m = img_recon_4d(:,:,ind_slice,ind_chunk,ind_dw) .* mask(:,:,ind_slice,ind_chunk,ind_dw);
                else                    
                    img_recon_m = img_recon_4d(:,:,ind_slice,ind_dw) .* mask(:,:,ind_slice);
                end
                
                if flag_resMethod==1
                    
                    %---------- Method 1 ----------
                    % Calculate final recon image size.
                    if strcmpi(DWIparams.fold_over_dir,'AP')
                        dwi_fov_p = DWIparams.FOV.AP;
                        if strcmpi(DWIparams.slice_orientation,'transverse')
                            dwi_fov_m = DWIparams.FOV.RL;
                        elseif strcmpi(DWIparams.slice_orientation,'sagittal')
                            dwi_fov_m = DWIparams.FOV.FH;
                        else
                            error('s_resize_DWI:main','Unknown DWI slice orientation')
                        end
                    elseif strcmpi(DWIparams.fold_over_dir,'RL')
                        dwi_fov_p = DWIparams.FOV.RL;
                        if strcmpi(DWIparams.slice_orientation,'transverse')
                            dwi_fov_m = DWIparams.FOV.AP;
                        elseif strcmpi(DWIparams.slice_orientation,'coronal')
                            dwi_fov_m = DWIparams.FOV.FH;
                        else
                            error('s_resize_DWI:main','Unknown DWI slice orientation')
                        end
                    else
                        error('s_resize_DWI:main','Unknown DWI fold-over directions')
                    end
                    
                    
                    % P_ORI.
                    % (epi_factor * shot * sense_factor) *
                    % RECONparams.ACQREC.ACQ_voxel_MPS(P_ORI) == N_resize_P_ORI *
                    % RECONparams.ACQREC.REC_voxel_MPS(P_ORI) But
                    % RECONparams.ACQREC.REC_voxel_MPS(P_ORI) is not adjusted here when
                    % scan user changed the recon voxel size to be isotropic as M_ORI.
                    % Then just use RECONparams.ACQREC.REC_voxel_MPS(M_ORI).
                    N_resize_P_ORI = round(ny * RECONparams.ACQREC.ACQ_voxel_MPS(2) / ...
                        RECONparams.ACQREC.REC_voxel_MPS(1));
                    N_recon_P_ORI = round(dwi_fov_p / RECONparams.ACQREC.REC_voxel_MPS(1));
                    
                    
                    % M_ORI.
                    N_resize_M_ORI = round(nx * RECONparams.ACQREC.ACQ_voxel_MPS(1) / ...
                        RECONparams.ACQREC.REC_voxel_MPS(1));
                    N_recon_M_ORI = round(dwi_fov_m / RECONparams.ACQREC.REC_voxel_MPS(2));
                    
                    
                    % Resize to recon voxel size.
                    % There is stripe-like artifact when complex image data is
                    % interpolated due to phase discontinuity. The artifact appears
                    % also when real and imaginary image data are interpolated
                    % separately. Then interpolate magnitude and phase data separately.
                    % This generate no artifact.
                    %                 img_recon_res_abs_m = imresize(abs(img_recon_m),[N_resize_P_ORI,size(img_recon_m,2)]);
                    %                 img_recon_res_phase_m = imresize(angle(img_recon_m),[N_resize_P_ORI,size(img_recon_m,2)]);
                    %                 img_recon_res_m = img_recon_res_abs_m.*exp(1i*img_recon_res_phase_m);
                    img_recon_res_abs_m = imresize(abs(img_recon_m),[N_resize_P_ORI,N_resize_M_ORI]);
                    img_recon_res_phase_m = imresize(angle(img_recon_m),[N_resize_P_ORI,N_resize_M_ORI]);
                    img_recon_res_m = img_recon_res_abs_m.*exp(1i*img_recon_res_phase_m);
                    
                else
                    
                    %---------- Method 2 ----------
                    % Calculate final recon image size.
                    
                    
                    % P_ORI.
                    % This uses the relationship as following,
                    % [ovs_res / ovs_fac * sense_fac * interp_fac = recon_res]
                    % and [FOV = recon_res * rec_vox_p] and also
                    % [acq_vox_p = ovs_fac * FOV / (epi_fac * shot * sense_fac)]
                    % then [interp_fac = acq_vox_p / rec_vox_p].
                    % Then the final resize matrix size can be
                    % [epi_fac * shot * sense_fac) / ovs_fac * interp_fac].
                    
                    %N_resize_P_ORI = ny/RECONparams.ENCima.oversample_factors(2)* ...
                    %    RECONparams.ENCima.interp_factors(2);
                    N_resize_P_ORI = floor(ny * RECONparams.ENCima.interp_factors(2));
                    N_recon_P_ORI = RECONparams.ENCima.recon_resolutions(2); % isotropic voxel
                    
                    
                    % M_ORI.
                    %N_resize_M_ORI = nx/RECONparams.ENCima.oversample_factors(1)* ...
                    %    RECONparams.ENCima.interp_factors(1);
                    N_resize_M_ORI = floor(nx * RECONparams.ENCima.interp_factors(1));
                    N_recon_M_ORI = RECONparams.ENCima.recon_resolutions(1);
                    
                    
                    % Resize to recon voxel size.
                    %                 img_recon_res_abs_m = imresize(abs(img_recon_m),[N_resize_P_ORI,size(img_recon_m,2)]);
                    %                 img_recon_res_phase_m = imresize(angle(img_recon_m),[N_resize_P_ORI,size(img_recon_m,2)]);
                    %                 img_recon_res_m = img_recon_res_abs_m.*exp(1i*img_recon_res_phase_m);
                    img_recon_res_abs_m = imresize(abs(img_recon_m),[N_resize_P_ORI,N_resize_M_ORI]);
                    img_recon_res_phase_m = imresize(angle(img_recon_m),[N_resize_P_ORI,N_resize_M_ORI]);
                    img_recon_res_m = img_recon_res_abs_m.*exp(1i*img_recon_res_phase_m);
                end
                
                
                % Determine zeropad or cut in resized recon image.
                if N_resize_P_ORI <= N_recon_P_ORI
                    
                    % P_ORI.
                    zeropad_vox_P_ORI = (N_recon_P_ORI - N_resize_P_ORI)/2;
                    zeropad_vox_pre_P_ORI = floor(zeropad_vox_P_ORI);
                    zeropad_vox_post_P_ORI = ceil(zeropad_vox_P_ORI);
                    
                    img_recon_res_zeropad_m = padarray(padarray( ...
                        img_recon_res_m,[zeropad_vox_pre_P_ORI,0],0,'pre'), ...
                        [zeropad_vox_post_P_ORI,0],0,'post');
                    
                    % M_ORI.
                    % Acquisition matrix size along M_ORI is usually larger
                    % than recon matrix size along M_ORI due to
                    % oversampling_factor of 2 along M_ORI. But sometimes recon
                    % matrix size is bigger than acquisition matrix size, when
                    % FOV is very small. So consider the case here along M_ORI.
                    %N_recon_M_ORI = N_recon_P_ORI; % due to isotropic recon image
                    if size(img_recon_res_m,2) >= N_recon_M_ORI
                        cut_vox_M_ORI = (size(img_recon_res_m,2) - N_recon_M_ORI)/2;
                        cut_vox_pre_M_ORI = floor(cut_vox_M_ORI);
                        cut_vox_post_M_ORI = ceil(cut_vox_M_ORI);
                        
                        img_recon_res_zeropad_cut_m = ...
                            img_recon_res_zeropad_m(:,cut_vox_pre_M_ORI+1:cut_vox_post_M_ORI+N_recon_M_ORI);
                        m = img_recon_res_zeropad_cut_m;
                    else
                        ovs_fac_M_ORI = RECONparams.ENCima.oversample_factors(1);
                        img_recon_res_zeropad_res_abs_m = imresize(abs(img_recon_res_zeropad_m), ...
                            [size(img_recon_res_zeropad_m,1),size(img_recon_res_zeropad_m,2)*ovs_fac_M_ORI]);
                        img_recon_res_zeropad_res_phase_m = imresize(angle(img_recon_res_zeropad_m), ...
                            [size(img_recon_res_zeropad_m,1),size(img_recon_res_zeropad_m,2)*ovs_fac_M_ORI]);
                        img_recon_res_zeropad_res_m = img_recon_res_zeropad_res_abs_m .* ...
                            exp(1i*img_recon_res_zeropad_res_phase_m);
                        
                        cut_vox_M_ORI = (size(img_recon_res_zeropad_res_m,2) - N_recon_M_ORI)/2;
                        cut_vox_pre_M_ORI = floor(cut_vox_M_ORI);
                        cut_vox_post_M_ORI = ceil(cut_vox_M_ORI);
                        
                        m = img_recon_res_zeropad_res_m(: , ...
                            cut_vox_pre_M_ORI+1:cut_vox_pre_M_ORI+N_recon_M_ORI);
                    end
                else
                    % P_ORI.
                    cut_vox_P_ORI = (N_resize_P_ORI - N_recon_P_ORI)/2;
                    cut_vox_pre_P_ORI = floor(cut_vox_P_ORI);
                    cut_vox_post_P_ORI = ceil(cut_vox_P_ORI);
                    
                    img_recon_res_cut_m = ...
                        img_recon_res_m(cut_vox_pre_P_ORI+1:cut_vox_pre_P_ORI+N_recon_P_ORI,:);
                    
                    % M_ORI.
                    %N_recon_M_ORI = N_recon_P_ORI; % due to isotropic recon image
                    if size(img_recon_res_m,2) >= N_recon_M_ORI
                        cut_vox_M_ORI = (size(img_recon_res_m,2) - N_recon_M_ORI)/2;
                        cut_vox_pre_M_ORI = floor(cut_vox_M_ORI);
                        cut_vox_post_M_ORI = ceil(cut_vox_M_ORI);
                        
                        img_recon_res_cut_cut_m = ...
                            img_recon_res_cut_m(:,cut_vox_pre_M_ORI+1:cut_vox_pre_M_ORI+N_recon_M_ORI);
                        m = img_recon_res_cut_cut_m;
                    else
                        ovs_fac_M_ORI = RECONparams.ENCima.oversample_factors(1);
                        img_recon_res_cut_res_abs_m = imresize(abs(img_recon_res_cut_m), ...
                            [size(img_recon_res_cut_m,1),size(img_recon_res_cut_m,2)*ovs_fac_M_ORI]);
                        img_recon_res_cut_res_phase_m = imresize(angle(img_recon_res_cut_m), ...
                            [size(img_recon_res_cut_m,1),size(img_recon_res_cut_m,2)*ovs_fac_M_ORI]);
                        img_recon_res_cut_res_m = img_recon_res_cut_res_abs_m .* ...
                            exp(1i*img_recon_res_cut_res_phase_m);
                        
                        cut_vox_M_ORI = (size(img_recon_res_cut_res_m,2) - N_recon_M_ORI)/2;
                        cut_vox_pre_M_ORI = floor(cut_vox_M_ORI);
                        cut_vox_post_M_ORI = ceil(cut_vox_M_ORI);
                        
                        m = img_recon_res_cut_res_m(:, ...
                            cut_vox_M_ORI+1:cut_vox_M_ORI+N_recon_M_ORI);
                    end
                end
                clear  img_recon_res_m  img_recon_res_zeropad_m  img_recon_res_zeropad_cut_m
                                
                
                %---------------------------------------------------------
                % Re-orient recon image data.
                
                %* Transverse: A-L-P-R in N-E-S-W
                %* Coronal: R-H-L-F in N-E-S-W
                %* Sagittal: A-F-P-H in N-E-S-W
                
                % Transverse.
                if strcmpi(DWIparams.slice_orientation,'transverse')
                    
                    if strcmpi(DWIparams.fold_over_dir,'AP')
                        
                        if strcmpi(DWIparams.fat_shift_dir,'P')
                            if strcmpi(DWIparams.scan_mode,'3D')
                                dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = m;
                            else                                
                                dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = m;
                            end
                            
                        elseif strcmpi(DWIparams.fat_shift_dir,'A')
                            m = fliplr(flipud(m));
                            if strcmpi(DWIparams.scan_mode,'3D')
                                dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = m;
                            else
                                dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = m;
                            end
                            
                        else
                            error('s_resize_DWI:main','Unknown ''DWIparams.fat_shift_dir''.')
                        end
                        
                    elseif strcmpi(DWIparams.fold_over_dir,'RL')
                        
                        if strcmpi(DWIparams.fat_shift_dir,'L')
                            m = flipud(permute(m,[2,1]));
                            if strcmpi(DWIparams.scan_mode,'3D')
                                dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = m;
                            else
                                dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = m;
                            end
                            
                        elseif strcmpi(DWIparams.fat_shift_dir,'R')
                            m = fliplr(permute(m,[2,1]));
                            if strcmpi(DWIparams.scan_mode,'3D')
                                dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = m;
                            else
                                dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = m;
                            end
                            
                        else
                            error('s_resize_DWI:main','Unknown ''DWIparams.fat_shift_dir''.')
                        end
                    else
                        error('s_resize_DWI:main','Unknown ''DWIparams.fold_over_dir''.')
                    end
                end
                
                
                % Sagittal.
                if strcmpi(DWIparams.slice_orientation,'sagittal')
                    
                    if strcmpi(DWIparams.fold_over_dir,'AP')
                        
                        if strcmpi(DWIparams.fat_shift_dir,'P')
                            % Check which one is correct.
                            m = fliplr(m);
                            %m = fliplr(flipud(m));
                            if strcmpi(DWIparams.scan_mode,'3D')
                                dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = m;
                            else
                                dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = m;
                            end
                            
                        elseif strcmpi(DWIparams.fat_shift_dir,'A')
                            m = fliplr(flipud(m));
                            if strcmpi(DWIparams.scan_mode,'3D')
                                dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = m;
                            else
                                dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = m;
                            end
                            
                        else
                            error('s_resize_DWI:main','Unknown ''DWIparams.fat_shift_dir''.')
                        end
                        
                    elseif strcmpi(DWIparams.fold_over,dir,'FH')
                        
                        if strcmpi(DWIparams.fat_shift_dir,'F')
                            warning('s_resize_DWI:main','SAGITTAL FH-F is not setup yet')
                            
                        elseif strmcpi(DWIparams.fat_shift_dir,'H')
                            warning('s_resize_DWI:main','SAGITTAL FH-H is not setup yet')
                            
                        else
                            error('s_resize_DWI:main','Unknown ''DWIparams.fat_shift_dir''.')
                        end
                    else
                        error('s_resize_DWI:main','Unknown ''DWIparams.fold_over_dir''.')
                    end
                end
                
                
                % Coronal.
                if strcmpi(DWIparams.slice_orientation,'coronal')
                    
                    if strcmpi(DWIparams.fold_over_dir,'RL')
                        
                        if strcmpi(DWIparams.fat_shift_dir,'L')
                            if strcmpi(DWIparams.scan_mode,'3D')
                                dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = m;
                            else
                                dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = m;
                            end
                            
                        elseif strcmpi(DWIparams.fat_shift_dir,'R')
                            warning('s_resize_DWI:main','CORONAL RL-R is not setup yet')
                            
                        else
                            error('s_resize_DWI:main','Unknown ''DWIparams.fat_shift_dir''.')
                        end
                        
                    elseif strcmpi(DWIparams.fold_over,dir,'FH')
                        
                        if strcmpi(DWIparams.fat_shift_dir,'F')
                            warning('s_resize_DWI:main','CORONAL FH-F is not setup yet')
                            
                        elseif strmcpi(DWIparams.fat_shift_dir,'H')
                            warning('s_resize_DWI:main','CORONAL FH-H is not setup yet')
                            
                        else
                            error('s_resize_DWI:main','Unknown ''DWIparams.fat_shift_dir''.')
                        end
                    else
                        error('s_resize_DWI:main','Unknown ''DWIparams.fold_over_dir''.')
                    end
                end
                
                
                
                %---------------------------------------------------------
                % Re-orient, again, into <default display orientation> before
                % tensor calculation using DWIparams.DW_GRAD which is in CRP
                % coordinate. See [f_get_dw_orientation.m] and E4p44.
                
                %* Transverse: A-L-P-R -> A-L-P-R
                %* Coronal: R-H-L-F -> H-L-F-R
                %* Sagittal: A-F-P-H -> H-P-F-A
                
                if strcmpi(DWIparams.slice_orientation,'transverse')
                    if strcmpi(DWIparams.scan_mode,'3D')
                        dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = m;
                    else
                        dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = m;
                    end
                elseif strcmpi(DWIparams.slice_orientation,'coronal')
                    if strcmpi(DWIparams.scan_mode,'3D')
                        dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = flipud(permute(m,[2,1]));
                    else
                        dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = flipud(permute(m,[2,1]));
                    end
                elseif strcmpi(DWIparams.slice_orientation,'sagittal')
                    if strcmpi(DWIparams.scan_mode,'3D')
                        dw_data_5d(:,:,ind_slice,ind_chunk,ind_dw,dw_ori+1) = permute(m,[2,1]);
                    else
                        dw_data_5d(:,:,ind_slice,ind_dw,dw_ori+1) = permute(m,[2,1]);
                    end
                else
                    error('s_resize_DWI:main','Unknown slice_orientation')
                end                
                
                % Clear.
                %clear  img_recon_m  img_recon_res*  m
                
            end % if RECONFLAGparams.res_DWI
            
        end % for ind_slice
        end % for ind_chunk
    end % for ind_dw
end % for dw_ori


% Save resized DWI images.
cd(saveReconDir_s)
if RECONFLAGparams.phaseCorr==0
    dw_data_no_corr_5d = dw_data_5d;
    %save  dw_data_no_corr_5d  dw_data_no_corr_5d
    if isempty(flag_average_idx)
        save  dw_data_no_corr_5d  dw_data_no_corr_5d
    else
        eval(sprintf('dw_data_no_corr_5d%s = dw_data_5d;',tag_average_idx))
        eval(sprintf('save  dw_data_no_corr_5d%s  dw_data_no_corr_5d%s',...
            tag_average_idx,tag_average_idx))
    end
end
if RECONFLAGparams.phaseCorr==1
    if RECONFLAGparams.resize_DWI==1
        if RECONFLAGparams.navFmapCorr==0 && RECONFLAGparams.navDeform==0
            save  dw_data_5d  dw_data_5d
        end
        if RECONFLAGparams.navFmapCorr==1
            dw_data_fmap_5d = dw_data_5d;
            %save  dw_data_fmap_5d  dw_data_fmap_5d
            if isempty(flag_average_idx)
                save  dw_data_fmap_5d  dw_data_fmap_5d
            else
                eval(sprintf('dw_data_fmap_5d%s = dw_data_5d;',tag_average_idx))
                eval(sprintf('save  dw_data_fmap_5d%s  dw_data_fmap_5d%s',...
                    tag_average_idx,tag_average_idx))
            end
        end
        if RECONFLAGparams.navDeform==1
            dw_data_deform_5d = dw_data_5d;
            save  dw_data_deform_5d  dw_data_deform_5d
        end
    else
        if RECONFLAGparams.navFmapCorr==0 && RECONFLAGparams.navDeform==0
            dw_data_no_res_5d = dw_data_5d;
            save  dw_data_no_res_5d  dw_data_no_res_5d
        end
        if RECONFLAGparams.navFmapCorr==1
            dw_data_fmap_no_res_5d = dw_data_5d;
            save  dw_data_fmap_no_res_5d  dw_data_fmap_no_res_5d            
        end
        if RECONFLAGparams.navDeform==1
            dw_data_deform_no_res_5d = dw_data_5d;
            save  dw_data_deform_no_res_5d  dw_data_deform_no_res_5d
        end
    end
end
clear  dw_data_5d  img_*  mask

% Pack memory.
pack

% end % if flag_resizeDWI

fprintf('\n\n')

% Show images.
% imageviewxd(dw_data_5d,'DW data')










