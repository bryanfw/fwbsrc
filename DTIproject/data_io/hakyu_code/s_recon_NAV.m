
%[s_recon_NAV] reconstructs NAV data for b-value, slices and shots.
%
%
% Last modified
% 2010.02.11.
% 2010.02.17. saveReconDir_s -> saveDataDir_s to separate DW data and recon data.
% 2010.02.18. Add code for kx and ky Hamming windowing before NAV recon.
% 2010.02.19. Add DATAFLAGparams.navWindowedRecon for windowing k-space NAV data 
%             before NAV recon.
% 2010.03.26. Try Hamming windowed recon along X as well as Y. Previously Y only
%             window has been used.
% 2010.07.12. Consider the case when foldover is 'RL'.
%             Put '% Match foldover direction.' part.
% 2010.07.22. Consider the case when foldover is 'RL'.
%             Put '% Adjust reconstructed object position variable with foldover and fatshift.' part.
% 2010.07.26. Change psi to PSI to remove a conflict to Matlab function
%             'psi()'.
% 2010.08.02. Don't use '% Adjust reconstructed object position variable
%             with foldover and fatshift.' part after recon. This will
%             require another adjustment in [Nav_fmapCorr.m] is run.
% 2010.08.09.
%   This is generated from [NAV_recon.m].
%   Don't use DATAFLAGparams.navWindowedRecon in future, because NAV recon
%   will be done w/ and w/o window as default.
%   Reconstruct with and without k-space window as default.
% 2010.08.17.
%   Add '%----- Get ky_range_v: method 2 -----'. This looks right than the
%   method 1. To be more accurate, NAV must be gridded onto ky position,
%   but sinc interpolation generates artifact after NAV recon. So gridding
%   is currently on hold.
% 2010.08.18.
%   Pack memory.
% 2010.09.21.
%   Use K_alias_6d and K_nav_6d instead of K_alias_5d and K_nav_5d
%   considering multiple non-zero b-values.
% 2010.09.22.
%   Consider REFparams.fold_over_dir='RL'.
% 2010.09.24.
%   Multiply (-1) to odd line of k-space (undersampled) data. This is only
%   the case for 3T.
% 2010.09.28.
%   Take care old format of k-space data, K_alias_5d and K_nav_5d.
% 2010.10.05.
%   Take input of GENparams to [f_sense_recon.m].
% 2010.12.07.
%   Multiply (-1) for 3T data.
% 2011.03.03.
%   Apply  'K_alias_shot_m(1:2:end,:) = K_alias_shot_m(1:2:end,:)*(-1);'
%   for 7T data when ZOOM is used.
% 2011.04.11.
%   Modifiy kx_cut = (nnx - nnx_new)/2; to kx_cut = round((nnx -
%   nnx_new)/2;) not to have errors.
% 2012.05.17.
%   Use [f_sense_recon_v2.m] for matching object position between SMAP and
%   NAV.
% 2012.05.23.
%   Take care of Volume coil case (REFparams.filename is empty).
% 2012.06.16.
%   Take care of scan_mode = 3D case.
% 2012.07.18.
%   gfactor is not saved to reduce memory and hdd space.
%
% Ha-Kyu



%% Recon NAV data for all DW and SLICES
warning off

% Load smap data.
if ~isempty(REFparams.filename)
    cd(saveDataDir_s)
    
    % Load some params struct. This may or may not include the necessary
    % fields based on how many times this script run.
    if ~isfield('DATAFLAGparams','matchObjPosDelta')
        cd(saveReconDir_s)
        if exist('DATAFLAGparams.matchObjPosDelta.mat','file')
            load('DATAFLAGparams.matchObjPosDelta.mat')
            DATAFLAGparams.matchObjPosDelta = matchObjPosDelta;
        end
    end
    if ~isfield('SMAPparams','filesize__I_smap_img_4d')
        cd(saveReconDir_s)
        if exist('SMAPparams.filesize__I_smap_img_4d.mat','file')
            load('SMAPparams.filesize__I_smap_img_4d.mat')
            SMAPparams.filesize__I_smap_img_4d = filesize__I_smap_img_4d;
        end
    end
    if ~isfield('SMAPparams','filesize__ref')
        cd(saveReconDir_s)
        if exist('SMAPparams.filesize__ref.mat','file')
            load('SMAPparams.filesize__ref.mat')
            SMAPparams.filesize__ref = filesize__ref;
        end
    end
    cd(saveDataDir_s)
    
    
    % Check if I_smap_img_4d is saved.
    % If I_smap_img_4d size is > 2GB, load this later.
    if exist('I_smap_img_4d.mat','file')
        load  I_smap_img_4d
        I_smap_nav_4d = I_smap_img_4d;  clear  I_smap_img_4d;
        load  I_mask_S_img_3d
        I_mask_S_nav_3d = I_mask_S_img_3d; clear  I_mask_S_img_3d;
    else
        % Load chunk-by-chunk SMAP and mask data below.
    end
    %load  I_smap_img_4d    
    %load  I_mask_S_nav_3d        
    cd(sharedDir_s)
    load  PSI
else
    PSI = 1; % Volume Coil case
end

% Match foldover direction.
if exist('I_smap_nav_4d','var')
    s_match_foldover_direction_NAV
end


% Load DATAFLAGparams.
%* This is when DATAFLAGparams.matchObjPosSMAP_IMG==1 and there is updated
%* DATAFLAGparams.matchObjPosDelta filed calculated in [s_recon_IMG3.m].
% if DATAFLAGparams.matchObjPosSMAP_IMG==1
%     cd(saveReconDir_s)
%     load('params','DATAFLAGparams')
%     fprintf('*** DATAFLAGparams is loaded ***\n')
% end



%% Recon NAV   
for dw_ori = 0:DWIparams.nDW_GRAD
    fprintf('\n\n')
    fprintf('RECON NAV DATA windowed or not, ori = [%d]\n',dw_ori)
        
    % Load data.
    cd(saveDataDir_s)
    eval(sprintf('load  dw_data_ori%.2d__R%dS%d',dw_ori,R,Ns))
        
    % Get data size.
    % Take care of old format, K_alias_5d.
    if exist('K_alias_6d','var')
        if strcmpi(DWIparams.scan_mode,'3D')
            [ny,nx,nz,nchunk,nc,navg,ndw] = size(K_alias_6d);
            [nny,nnx,nnz,nnchunk,nnc,nnavg,nndw] = size(K_nav_6d);
            clear  K_alias_6d
            % Report.
            fprintf('    K_nav_6d,                    [nny,nnx,nnz,nnchunk,nnc,nnavg,nndw]      = [%d,%d,%d,%d,%d,%d,%d].\n', ...
                nny,nnx,nnz,nnchunk,nnc,nnavg,nndw)
        else
            [ny,nx,nz,nc,navg,ndw] = size(K_alias_6d);
            [nny,nnx,nnz,nnc,nnavg,nndw] = size(K_nav_6d);
            nchunk = 1;
            nnchunk = 1;
            clear  K_alias_6d
            % Report.
            fprintf('    K_nav_6d,                    [nny,nnx,nnz,nnc,nnavg,nndw]      = [%d,%d,%d,%d,%d,%d].\n', ...
                nny,nnx,nnz,nnc,nnavg,nndw)
        end
    else
        [ny,nx,nz,nc,navg] = size(K_alias_5d);
        [nny,nnx,nnz,nnc,nnavg] = size(K_nav_5d);
        ndw = 1;
        nndw = 1;
        K_alias_6d = reshape(K_alias_5d,[ny,nx,nz,nc,navg,ndw]);
        K_nav_6d = reshape(K_nav_5d,[nny,nnx,nnz,nnc,nnavg,nndw]);
        clear  K_alias_5d  K_nav_5d
        
        [ny,nx,nz,nc,navg,ndw] = size(K_alias_6d);
        [nny,nnx,nnz,nnc,nnavg,nndw] = size(K_nav_6d);
        clear  K_alias_6d
        
        % Report.
        fprintf('    K_nav_6d,                    [nny,nnx,nnz,nnc,nnavg,nndw]      = [%d,%d,%d,%d,%d,%d].\n', ...
            nny,nnx,nnz,nnc,nnavg,nndw)        
    end
    pack
    fprintf('    Memory packed\n')
    
    
    %% Prepare NAV data
    %*** This comes from [prepare_IMG_NAV__scan********_r****.m]
    
    % Calculate NAV k-space data position in IMG k-space data matrix.
    
    %----- Get ky_range_v: method 1 -----
%     ky_img_v = -floor(ny/2):(ceil(ny/2)-1);
%     ky_nav_v = -floor(DWIparams.EPI_FACTOR(2)/2) : (ceil(DWIparams.EPI_FACTOR(2)/2) - 1);
%     ky_nav_begin = find(ky_img_v==ky_nav_v(1)); % counted from bottom of k-space image (not index)
%     ky_nav_end = find(ky_img_v==ky_nav_v(end));
%     ky_nav_bot = ny - ky_nav_begin + 1; % bottom index of starting nav acq
%     ky_nav_top = ny - ky_nav_end + 1;   % top index of ending nav acq
%     ky_range_v = ky_nav_bot:-1:ky_nav_top; % range of ky lines to
    
    %----- Get ky_range_v: method 2 -----
    % Another: index of ky_nav.
    ky_img_v = -floor(ny/2):(ceil(ny/2)-1);
    ky_nav_v = -floor(DWIparams.EPI_FACTOR(2)/2) : (ceil(DWIparams.EPI_FACTOR(2)/2) - 1);
    ind_zero_img = find(ky_img_v==0);
    ind_zero_nav = find(ky_nav_v==0);
    ky_nav_top = ind_zero_img - floor(DWIparams.EPI_FACTOR(2)/2);
    ky_nav_bot = ind_zero_img + (ceil(DWIparams.EPI_FACTOR(2)/2)-1);
    ky_range_v = ky_nav_bot:-1:ky_nav_top;
        
    
    % Apply k-space window to reduce Gibb's ringing.    
    %* Generate Hamming window to reduce Gibb's ringing.
    %* Use both of ky and kx (after cutting off to generate isotropic
    %* in-plane resolution) direction.
    
    %* Match the ratio of voxel size in x and y in NAV as is in IMG by
    %* calculating f, an x directional voxel size of NAV matching
    %* acq_vox_p(IMG)/acq_vox_m(IMG) = acq_vox_p(NAV)/acq_vox_m(NAV)
    %* f = RECONparams.nav__acq_vox_p * RECONparams.acq_vox_m / ...
    %*   RECONparams.acq_vox_p; % new x voxel size to be
    
    nav__acq_vox_p = RECONparams.ACQREC.ACQ_voxel_MPS(2) * ...
        RECONparams.ACQREC.number_of_shots;
    f = nav__acq_vox_p * RECONparams.ACQREC.ACQ_voxel_MPS(1) / ...
        RECONparams.ACQREC.ACQ_voxel_MPS(2); % f must be the number of shots.
    
    nnx_new = nnx/f; % new x matrix size to be
    if mod(nnx-floor(nnx_new),2)==0
        nnx_new = floor(nnx_new);
    else
        nnx_new = ceil(nnx_new);
    end
    kx_cut = round((nnx - nnx_new)/2);
    kx_range_v = kx_cut+1:kx_cut+nnx_new;
    
    hwin_m = zeros(ny,nnx,'single');    
    %hx_v = hamming(nnx)';   % don't match the voxel size ratio
    hx_v = hamming(length(kx_range_v))';   % match the voxel size ratio    
    hy_v = hamming(length(ky_nav_v));
    
    hwin = hy_v*hx_v;    
    %hwin_m(ky_range_v,:) = hwin;    % don't match the voxel size ratio
    hwin_m(ky_range_v,kx_range_v) = hwin;  % match the voxel size ratio
    
    
    %% Read and reconstruct data
    if dw_ori==0
        win_v = 0:1;
        %win_v = 1;
    else
        win_v = 1;
    end       
    for ind_windowed = win_v
        if strcmpi(DWIparams.scan_mode,'3D')
            %I_nav_sense_7d = zeros(ny,nnx,nnz,nnchunk,nnc,nnavg,nndw,Ns); % nav into img size
            I_nav_sense_7d = zeros(ny,nnx,nnz,1,nnc,nnavg,nndw,Ns); % keep chunk dim to 1
        else
            I_nav_sense_7d = zeros(ny,nnx,nnz,nnc,nnavg,nndw,Ns); % nav into img size
        end
        for ind_chunk = 1:nnchunk
            if strcmpi(DWIparams.scan_mode,'3D')
                I_nav_sense_7d = I_nav_sense_7d*0;
            end
            for ind_coil = 1:nnc
                for ind_sl = 1:nnz
                    for ind_avg = 1:nnavg
                        for ind_dw = 1:nndw
                            if strcmpi(DWIparams.scan_mode,'3D')
                                K_nav_m = K_nav_6d(:,:,ind_sl,ind_chunk,ind_coil,ind_avg,ind_dw);
                            else
                                K_nav_m = K_nav_6d(:,:,ind_sl,ind_coil,ind_avg,ind_dw);
                            end
                            for ind_shot = 1:Ns
                                % Take each shot nav data. From bottom to top of k-space image.
                                ky_shot_v = nny+1-ind_shot:-Ns:1;
                                
                                % Apply window or not.
                                K_nav_shot_m = zeros(ny,nnx);
                                if ind_windowed == 0
                                    K_nav_shot_m(ky_range_v,:) = ...
                                        K_nav_m(ky_shot_v,:);
                                else
                                    K_nav_shot_m(ky_range_v,kx_range_v) = ...
                                        K_nav_m(ky_shot_v,kx_range_v);
                                    K_nav_shot_m = K_nav_shot_m .* hwin_m;
                                end
                                
                                %---------- START ----------
                                % Need this for 3T.
                                if GENparams.B0==30000
                                    if ~strcmpi(GENparams.coilID,'SENSE-Breast-4')
                                        % [scan20100921_r2111__Smith_4117_3T] must multiply (-1).
                                        K_nav_shot_m(1:2:end,:) = K_nav_shot_m(1:2:end,:)*(-1); % original
                                        
                                        % [scan20101026_r2359__Smith_2402_3T] shouldn't multiply.
                                        % [scan20101207_r2359__Smith_4202_3T] should apply
                                    end
                                elseif GENparams.B0==70000
                                    if DWIparams.ZOOM==true
                                        K_nav_shot_m(1:2:end,:) = K_nav_shot_m(1:2:end,:)*(-1); % original
                                    end
                                    if strcmpi(GENparams.coilID,'RX-Intf-1_Quad-TR-1')
                                        K_nav_shot_m(1:2:end,:) = K_nav_shot_m(1:2:end,:)*(-1);
                                        % [scan20120420_r4140__Anderson_306871]
                                        % [scan20120515_r4140__Jeong_999999]
                                        % need *(-1)
                                    end
                                end
                                %---------- END ----------
                                
                                % Just use original size data.
                                I_nav_shot_m = ift2(K_nav_shot_m);
                                if strcmpi(DWIparams.scan_mode,'3D')
                                    %I_nav_sense_7d(:,:,ind_sl,ind_chunk, ...
                                    %    ind_coil,ind_avg,ind_dw,ind_shot) = ...
                                    %    I_nav_shot_m;
                                    I_nav_sense_7d(:,:,ind_sl,1, ...
                                        ind_coil,ind_avg,ind_dw,ind_shot) = ...
                                        I_nav_shot_m;
                                else
                                    I_nav_sense_7d(:,:,ind_sl,ind_coil,ind_avg,ind_dw,ind_shot) = ...
                                        I_nav_shot_m;
                                end
                            end % ind_shot
                        end % ind_dw
                    end % ind_avg
                end % for ind_sl
            end % for ind_coil
            
            % Save nav as chunk.
            cd(saveDataDir_s)
            if strcmpi(DWIparams.scan_mode,'3D')
                if ~exist(sprintf('I_nav_sense_chunk%.2d_7d.mat',ind_chunk),'file')
                    eval(sprintf('save  I_nav_sense_chunk%.2d_7d.mat  I_nav_sense_7d',ind_chunk))
                    fprintf('      I_nav_sense_chunk%.2d_7d.mat is saved\n',ind_chunk)
                else
                    fprintf('      I_nav_sense_chunk%.2d_7d.mat exists\n',ind_chunk)
                end
            end
            
        end % ind_chunk
        clear  K_nav_m  K_nav_shot_m  I_nav_shot_m  ky_shot_v
        pack
        
        % Get data size.
        if strcmpi(DWIparams.scan_mode,'3D')
            [snny,snnx,snnz,snnchunk,snnc,snnavg,snndw,snns] = size(I_nav_sense_7d);
            snnchunk = DWIparams.chunks;
            nav_size_v = [snny,snnx,snnz,snnchunk,snnc,snnavg,snndw,snns];
            
            % Report.
            fprintf('\n')
            fprintf('    I_nav_sense_7d, windowed[%d], [snny,snnx,snnz,snnchunk,snnc,snnavg,snndw,snns]  = [%d,%d,%d,%d,%d,%d,%d,%d]\n', ...
                ind_windowed,snny,snnx,snnz,snnchunk,snnc,snnavg,snndw,snns)
        else
            [snny,snnx,snnz,snnc,snnavg,snndw,snns] = size(I_nav_sense_7d);
            snnchunk = 1;
            
            % Report.
            fprintf('\n')
            fprintf('    I_nav_sense_7d, windowed[%d], [snny,snnx,snnz,snnc,snnavg,snndw,snns]  = [%d,%d,%d,%d,%d,%d,%d]\n', ...
                ind_windowed,snny,snnx,snnz,snnc,snnavg,snndw,snns)
        end        
        
        % Generate dummy I_smap_nav_4d and I_mask_S_nav_4d for Volume Coil
        % case.
        if isempty(REFparams.filename)
            I_smap_nav_4d = ones(snny,snnx,snnz*snnchunk,snnc);
            I_mask_S_nav_3d = ones(snny,snnx,snnz*snnchunk,snnc);
        end
        
        % Generate dummy I_smap_nav_4d when I_smap_nav_4d couldn't be saved
        % due to its filesize.
        if ~exist('I_smap_nav_4d','var')
            %I_smap_nav_4d =
            %zeros(snny*R,snnx,snnz*snnchunk,snnc,'single'); % Takes too
            %much memory
            
            I_smap_nav_4d = 0; % it will be taken care in [f_sense_recon_v2.m]
            I_mask_S_nav_3d = 0;
        end
        clear  snny  snnx  snnz  snnc  snnavg  snndw  snns
        
        % Generate dummy I_nav_sense_7d for 3D acquisition since each chunk
        % data is saved in above routine. Then input for this data will be
        % just scalar 0.
        if strcmpi(DWIparams.scan_mode,'3D')
            I_nav_sense_7d = 0;
        end
        
        % Recon all slices of NAV.
        %[nav_recon_6d,nav_gfactor_6d] = f_sense_recon(I_nav_sense_7d, I_smap_nav_4d, ...
        %    PSI, I_mask_S_nav_3d, R, SMAPparams, GENparams, dw_ori);
        
        % $$$ verify if SMAP slice position matches NAV kz position for each chunk $$$
        
        % For 3D acquisition, I_nav_sense_7d has single chunk. All chunk
        % data must be loaded inside the function below.
        
        if ind_windowed==0
            fname_s = sprintf('recon_nav_ori%.2d__R%dS%d',dw_ori,R,Ns);
        else
            fname_s = sprintf('recon_nav_win_ori%.2d__R%dS%d',dw_ori,R,Ns);
        end
        
        DATAFLAGparams.saveDataDir_s = saveDataDir_s; % this is temporary   
        
        cd(saveReconDir_s)
        if ~exist(sprintf('%s.mat',fname_s),'file')
            fprintf('      %s.mat doesn''t exist. Then recon\n',fname_s)
            cd(saveDataDir_s)
            if strcmpi(DWIparams.scan_mode,'3D')
                [nav_recon_6d,nav_gfactor_6d] = f_sense_recon_v3(I_nav_sense_7d, ...
                    I_smap_nav_4d, PSI, I_mask_S_nav_3d, R, SMAPparams, ...
                    GENparams, DATAFLAGparams, REFparams, DWIparams, dw_ori, nav_size_v);
            else
                [nav_recon_6d,nav_gfactor_6d] = f_sense_recon_v2(I_nav_sense_7d, ...
                    I_smap_nav_4d, PSI, I_mask_S_nav_3d, R, SMAPparams, ...
                    GENparams, DATAFLAGparams, REFparams, DWIparams, dw_ori);
            end
        else
            fprintf('      %s.mat exist. Then skip recon\n',fname_s)
        end
        
        
        % Save NAV data.
        cd(saveReconDir_s)
        if ~exist(sprintf('%s.mat',fname_s),'file')
            if ind_windowed==0
                fname_s = sprintf('recon_nav_ori%.2d__R%dS%d',dw_ori,R,Ns);
                %eval(sprintf('save  %s  nav_recon_6d  nav_gfactor_6d',fname_s))
                eval(sprintf('save  %s  nav_recon_6d',fname_s))
            else
                nav_recon_win_6d = nav_recon_6d;
                nav_gfactor_win_6d = nav_gfactor_6d;
                clear  nav_recon_6d  nav_gfactor_6d
                fname_s = sprintf('recon_nav_win_ori%.2d__R%dS%d',dw_ori,R,Ns);
                %eval(sprintf('save  %s  nav_recon_win_6d  nav_gfactor_win_6d',fname_s))
                eval(sprintf('save  %s  nav_recon_win_6d',fname_s))
            end
            fprintf('\n')
            fprintf('    [%s] is saved at \n',fname_s)
            fprintf('    [%s]\n',saveReconDir_s)
        else
            fprintf('\n')
            fprintf('    [%s] exists at \n',fname_s)
            fprintf('    [%s]\n',saveReconDir_s)
        end
        
        clear  nav_recon_6d  nav_recon_win_6d  nav_gfactor_6d
        clear  I_nav_sense_7d
        pack
        
        % Delete nav chunks.
        cd(saveDataDir_s)
        if strcmpi(DWIparams.scan_mode,'3D')
%             for ind_chunk = 1:nnchunk
%                 if exist(sprintf('I_nav_sense_chunk%.2d_7d.mat',ind_chunk),'file')
%                     eval(sprintf('delete  I_nav_sense_chunk%.2d_7d.mat',ind_chunk))
%                     fprintf('      I_nav_sense_chunk%.2d_7d.mat is deleted\n',ind_chunk)
%                 end
%             end
            delete  I_nav_sense_chunk*.mat
        end
        
    end % for ind_windowed
    
    % Clear.
    clear  K_alias*  K_nav*
    pack
    
end % for dw_ori
clear  fname_s  I_smap_nav_4d  I_mask_S_nav_3d  nav_recon_*  PSI

% Pack memory.
pack

fprintf('\n\n')



%% END





