
%[s_recon_DWI2] reconstructs IMG with NAV phase correction.
%
%
% See s_recon_DWI2_pct, test_s_recon_time
%
% Last modified
% 2012.01.20
%   Reconstruction routine is modified for faster reconstruction. By using
%   the routine, recon time reduced by half. This routine is from
%   [test_s_recon_DWI__precalc.m] and [test_s_precalc.m] which compare
%   algorithms between precalculation of encoding matrix and conventional
%   columnwise reconstruction both with improved algorithm for calculating
%   the encoding matrix.
%
%   It seems that there is no real advantage using precalculated encoding
%   matrix over columnwise reconstruction once the improved routine is
%   used. But the precalculation should be useful for MR scanner for 'recon
%   mode' of 'realtime' or 'immediate' is used since the computational
%   burden is much reduced for preparing encoding matrix before inversion
%   of the matrix.
%
% 2012.01.23
%   Recon time is compared,
%       with pre-calculation, [283.863+79.404]sec
%       without pre-calcution, with sub-encoding matrix, [329.339]sec
%       without pre-calculation, without sub-encoding matrix, [315.451]sec
%   Then use without pre-calculatio and without sub-encoding matrix.
%   However, as described on 2012.01.20, the use of pre-calculation must be
%   used for 'realtime' image recon.
%   For using this function, all necessary parameters are described in
%   [s_batch__scan20111122_r3314__Anderson_305832.m], and all other batch
%   script must be modified accrodingly.
%
%   See also s_batch__scan20111122_r3314__Anderson_305832, s_recon_DWI,
%   test_s_recon_time, test_s_recon_time_parallel
%
% 2012.03.10.
%   'flag_average_idx' is added to selectively recon for corresponding
%   average index. If it is [], then use all averages.
%
% 2012.05.23.
%   Take care of Volume Coil case (REFparams.filename is empty).
% 2012.06.13.
%   Take care of 3D acquisition for nz update. Don't use DWIparams.nSLICE
%   for number of slice because it doesn't show oversampled number of
%   slices along kz in 3D acquisition. Rather use
%   RECONparams.ENCima.oversample_resolutions(3) for number of slices along
%   kz.
% 2012.06.16.
%   Take care of 3D multi-chunk DTI case.
% 2012.07.05.
%   Keep Head-to-Foot slice order of SMAP by setting
%   SMAPparams.flip_slice_order == false for 3D case. Then regenerate
%   ind_slice1 accordingly.
% 2012.07.18.
%   g-factor is not generated to reduce memory and hdd space.
%
% Ha-Kyu



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
fprintf('RECONFLAGparams\n')
disp(RECONFLAGparams)



%% SENSE recon DWI for all b-values and slices

warning off

% Recon even when no phase correction is applied, because s_recon_IMG.m
% only recon images when b = 0. In here with no phase correction, all navg
% will be used for recon.

%if (RECONFLAGparams.phaseCorr)==1 && (BATCHparams.recon_data==1)
if BATCHparams.recon_data==1
    
    % Load SMAP for IMG recon.
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
            load  I_mask_S_img_3d
        else
            load  I_smap_4d
            [Ny,Nx,Nz,Nc] = size(I_smap_4d);
            %Nz = DWIparams.chunks * DWIparams.kz;
            Nz = DWIparams.kz; % [s_reslice_map_3d_v3.m]
            clear  I_smap_4d
        end
        %load  I_smap_img_4d        
        cd(sharedDir_s)
        load PSI
    else
        % Volume coil case.
        PSI = 1;
    end
    
    % Match foldover direction.
    % This is to match the object position of SMAP to foldover and fat shift
    % direction of DWI.
    if ~isempty(REFparams.filename)
        cd(saveDataDir_s)
        if exist('I_smap_img_4d','var')
            s_match_foldover_direction_IMG
            [Ny,Nx,Nz,Nc] = size(I_smap_img_4d);
        end        
    end
    
    % Take care of Volume Coil case (REFparams.filename is empty).
    if isempty(REFparams.filename)
        cd(saveReconDir_s)
        load(sprintf('recon_img_ori00__R%dS%d.mat',R,Ns))
        if strcmpi(DWIparams.scan_mode,'3D')
            [Ny,Nx,Nz,Nchunk,Navg] = size(img_recon_5d);
            clear  img_recon_5d
            Nc = DWIparams.nCOIL;
            I_smap_img_4d = ones(Ny,Nx,Nz*Nchunk,Nc);
            I_mask_S_img_3d = ones(Ny,Nx,Nz*Nchunk);
        else
            [Ny,Nx,Nz,Navg] = size(img_recon_5d);
            clear  img_recon_5d
            Nc = DWIparams.nCOIL;
            I_smap_img_4d = ones(Ny,Nx,Nz,Nc);
            I_mask_S_img_3d = ones(Ny,Nx,Nz);
        end
    end
    cd(saveDataDir_s)
    
    % Get image-space mixing matrix.
    A_c = cell(1,Ns);
    for ind_shot = (1:Ns)
        A_m = [];
        for ind_fold = 1:R*Ns
            A_m = [A_m, eye(Ny/(R*Ns))*exp(1i*2*pi*(ind_fold-1)*(ind_shot-1)/Ns)];
        end
        A_c{1,ind_shot} = sparse(A_m);
    end
    
    
    % Construct noise correlation matrix.
    %* This is done by the Cholesky factorization of the original noise
    %* correltation matrix to decorrelate coil sensitivity map and aliased
    %* image data.
    %* RECONFLAGparams.psiMethod of 1 and 2 generates the same results.
    
    %* Take care of volume coil case (REFparams.filename is empty).
    if isempty(REFparams.filename)
        RECONFLAGparams.psiMethod = 0;
        PSI = 1;
    end
    
    if RECONFLAGparams.psiMethod==1
        %*** Method 1: E_prime(j) = PHI * E(j,:,s)
        L = chol(PSI,'lower');
        PP = kron(inv(L),eye(Ny/(R*Ns)));
        clear  L  PSI
    elseif RECONFLAGparams.psiMethod==2
        %*** Method 2: E_prime(j) = PHI * E(j)
        L = chol(PSI,'lower');
        PP = kron(inv(L),eye(Ny/R));
        clear  L  PSI
    elseif RECONFLAGparams.psiMethod==0
        %*** Method 0: Do nothing
        PP = kron(eye(size(PSI)),eye(Ny/(R*Ns)));
        clear  PSI
    else
        error('s_recon_DWI:main','Unknown RECONFLAGparams.psiMethod.')
    end
    
    
    % Report RECONFLAGparams.
    %disp(RECONFLAGparams)
    
    
    % SENSE reconstruction of DWI IMG data with phase correction.
    for dw_ori = 0:DWIparams.nDW_GRAD
        fprintf('Recon ORI[%d], RECONFLAGparams.psiMethod[%d]\n', ...
            dw_ori,RECONFLAGparams.psiMethod)
        
        % Load I_alias_sense_7d, if exist.
        % This is the case when [I_alias_sense_7d.mat] is small and saved
        % as a whole-file, not a slice-file as modified on 2010.11.07.
        % See [s_recon_IMG.m].
        cd(saveDataDir_s)
        fname_img_s = sprintf('I_alias_sense_ori%.2d__R%dS%d',dw_ori,R,Ns);
        if exist([fname_img_s,'.mat'],'file')
            % This is a whole-file case.
            
            % Load the data.
            load(fname_img_s)
            
            % Take care of old format, I_alias_sense_6d.
            % I_alias_sense_6d is always a whole-file.
            if exist('I_alias_sense_6d','var')
                [ny,nx,nz,nc,navg,ns] = size(I_alias_sense_6d);
                ndw = 1; % in this case, number of b~=0 is always 1
                I_alias_sense_7d = reshape(I_alias_sense_6d,[ny,nx,nz,nc,navg,ndw,ns]);
                clear  I_alias_sense_6d
            end
            
            % Get the data size.
            [ny,nx,nz,nc,navg,ndw,ns] = size(I_alias_sense_7d);
            
            % Set flag.
            flag__I_alias_sense_7d__wholefile = 1;
        else
            % This is a slice-file case.
            if strcmpi(DWIparams.scan_mode,'3D')                
                
                fname_img_s = sprintf('I_alias_sense_sl01_chunk01_ori%.2d__R%dS%d', ...
                    dw_ori,R,Ns); % check for slice 01 chunk 01
                if exist(sprintf('%s.mat',fname_img_s),'file')                    
                
                    % Load the data.
                    load(fname_img_s)
                    
                    % Get the data size: nz and nchunk are always 1.
                    [ny,nx,nz,nchunk,nc,navg,ndw,ns] = size(I_alias_sense_sl01_chunk01_7d);                                        
                else
                    fname_img_s = sprintf('H_alias_sense_sl01_chunk01_ori%.2d__R%dS%d', ...
                        dw_ori,R,Ns); % check for slice 01 chunk 01
                    
                    % Load the data.
                    load(fname_img_s)                    
                    
                    % Get the data size: nz and nchunk are always 1.
                    [ny,nx,nz,nchunk,nc,navg,ndw,ns] = size(H_alias_sense_sl01_chunk01_7d);
                end
                
            else
                fname_img_s = sprintf('I_alias_sense_sl01_ori%.2d__R%dS%d', ...
                    dw_ori,R,Ns); % check for slice 01
                
                % Load the data.
                load(fname_img_s)
                
                % Get the data size.
                [ny,nx,nz,nc,navg,ndw,ns] = size(I_alias_sense_sl01_7d);
                nchunk = 1;
            end
            
            
            % Update slice number.
            %* For 3D acquisition, there is kz oversampling factor. Then
            %* don't update nz using DWIparams.nSLICE. DWIparams.nSLICE
            %* doesn't show oversampled number of slices (along kz) in 3D.
            %* Rather use RECONparams.ENCima.oversample_resolutions(3).
            if strcmpi(DWIparams.scan_mode,'3D')
                nz = RECONparams.ENCima.oversample_resolutions(3);
                nchunk = DWIparams.chunks;
            else
                nz = DWIparams.nSLICE;
                nchunk = DWIparams.chunks;
            end
            
            % Set flag.
            flag__I_alias_sense_7d__wholefile = 0;
            
            % Clear data.
            clear  I_alias_sense*  H_alias_sense*
        end
        
        
        % Load NAV recon data.
        cd(saveReconDir_s)
        fname_nav_s = sprintf('%s_ori%.2d__R%dS%d',fname_nav_head_s,dw_ori,R,Ns);
        
        
        %----- TEMP -----
        flag = 1;
        while flag
            if ~exist(sprintf('%s.mat',fname_nav_s),'file')
                fprintf('    wait for 10 minutues\n')
                pause(60*10)
                flag = 1;
            else
                s=dir(sprintf('%s.mat',fname_nav_s));
                if s.bytes > 200*10^6
                    fprintf('    do final recon\n')
                    flag = 0;
                else
                    flag=1;
                end
            end
        end
        %sendEmail(sprintf('%s is loaded for final recon for [%s]', ...
        %    fname_nav_s,DWIparams.filename),' ')
        %----- TEMP -----
            
            
        load(fname_nav_s)
        eval(sprintf('nav_recon_6d = %s;',dataname_nav_s));
        if strcmpi(DWIparams.scan_mode,'3D')
            [nny,nnx,nnz,nnchunk,nnavg,nndw,nns] = size(nav_recon_6d);
        else
            [nny,nnx,nnz,nnavg,nndw,nns] = size(nav_recon_6d);
        end
        if ~strcmpi(dataname_nav_s,'nav_recon_6d')
            eval(sprintf('clear  %s',dataname_nav_s))
        end
        clear  nav_gfactor*
        
        
        % Clear NaN.
        nav_recon_6d(isnan(nav_recon_6d))=0;
        
        
        % Adjust average index for image reconstruction.
        if isempty(flag_average_idx)
            avg_v = 1:navg; % recon using all navg images
        else
            avg_v = flag_average_idx; % recon for one of navg image
            navg = length(avg_v);
        end
        fprintf('  [average index,navg]\n'), disp([avg_v,navg])
        
        
        % Reserve DWI IMG recon output.
        %* Because all 'navg' are combined to one during recon, actual
        %* recon data is only 4D considering multiple nonzero b-values.
        if strcmpi(DWIparams.scan_mode,'3D')
            img_recon_4d = zeros(Ny,nx,nz,nchunk,ndw,'single');
            %img_gfactor_4d = zeros(Ny,nx,nz,nchunk,ndw,'single');
        else
            img_recon_4d = zeros(Ny,nx,nz,ndw,'single');
            %img_gfactor_4d = zeros(Ny,nx,nz,ndw,'single');
        end
        
        % TEMP-----
        %if dw_ori==3
        % TEMP-----
        
        % Get slice vector.
        slice_v = [];
        ref = floor(nz/2)+1;
        slice_v(1) = ref;
        for ind = 1:nz-1
            slice_v(ind+1) = slice_v(ind) + ind * (-1)^ind;
        end
        
        % Recon each slice.
        for ind_dw = 1:ndw
            for ind_chunk = 1:nchunk                
                for ind_slice = slice_v
                    
                    % nz is slice/chunk for 3D case. But SMAP and FMAP slices
                    % are 1:nkz*nchunk. Then make new index for SMAP slices.                    
                    if strcmpi(DWIparams.scan_mode,'3D')
                        %ind_slice1 = ind_slice + nz * (ind_chunk - 1);
                        %ind_slice1 = (nz - ind_slice + 1) + nz * (ind_chunk - 1);
                        
                        % For [s_reslice_smap_v2.m] with unreversed SMAP slice
                        % order prescribed in [s_get_SMAPparams.m] and
                        % [f_coilsensemap.m]. 2012.07.05.
                        %ind_slice1 = (nchunk-ind_chunk)*nz + ind_slice;
                        
                        % For [s_reslice_smap_v3.m].
                        ind_slice1 = ind_slice;
                    else
                        ind_slice1 = ind_slice;
                    end
                    
                                        
                    % Load I_smap_img_sl%.2d_4d for each slice.
                    % If I_smap_img_4d size is > 2GB, load this.
                    cd(saveDataDir_s)
                    if (~exist('I_smap_img_4d.mat','file'))                        
                        eval(sprintf('load  I_smap_img_chunk%.2d_4d',ind_chunk))
                        eval(sprintf('load  I_mask_S_img_chunk%.2d_3d',ind_chunk))
                        ind_slice2 = ind_slice1; % always 1
                        
                        s_match_foldover_direction_IMG
                    else
                        ind_slice2 = ind_slice1;
                    end
                    cd(saveReconDir_s)
                    
                    
                    fprintf('\n')
                    % Load I_alias* data or set slice_idx.
                    if flag__I_alias_sense_7d__wholefile==0
                        % Load I_alias* data when this is a slice-file case.
                        cd(saveDataDir_s)
                        if strcmpi(DWIparams.scan_mode,'3D')
                            
                            % All dimensions are in image-space.
%                             fname_img_s = sprintf('I_alias_sense_sl%.2d_chunk%.2d_ori%.2d__R%dS%d', ...
%                                 ind_slice,ind_chunk,dw_ori,R,Ns);
%                             load(fname_img_s)
%                             eval(sprintf('I_alias_sense_7d = I_alias_sense_sl%.2d_chunk%.2d_7d;', ...
%                                 ind_slice,ind_chunk))
                            
                            % Hybrid-space in in-plane dimensions. IFT
                            % along kz direction is required later.
                            fname_img_s = sprintf('H_alias_sense_sl%.2d_chunk%.2d_ori%.2d__R%dS%d', ...
                                ind_slice,ind_chunk,dw_ori,R,Ns);
                            load(fname_img_s)
                            eval(sprintf('I_alias_sense_7d = H_alias_sense_sl%.2d_chunk%.2d_7d;', ...
                                ind_slice,ind_chunk)) % keep the name I_alias_sense_7d
                            fprintf('  [%s] is loaded\n',fname_img_s)
                        else                            
                            fname_img_s = sprintf('I_alias_sense_sl%.2d_ori%.2d__R%dS%d', ...
                                ind_slice,dw_ori,R,Ns);
                            load(fname_img_s)
                            eval(sprintf('I_alias_sense_7d = I_alias_sense_sl%.2d_7d;', ...
                                ind_slice)) % keep the name I_alias_sense_7d
                            fprintf('  [%s] is loaded\n',fname_img_s)
                        end
                        %eval(sprintf('clear  I_alias_sense_sl%.2d_7d',ind_slice))
                        clear  I_alias_sense_sl*
                        cd(saveReconDir_s)
                        
                        % Set slice_idx for use instead of ind_slice. Because
                        % there is only one slice in this case.
                        % Add chunk_idx.
                        slice_idx = 1;
                        chunk_idx = 1;
                    else
                        % This is whole-file case.
                        
                        % Set slice_idx.
                        slice_idx = ind_slice;
                    end
                    
                    % Shift IMG data.
                    if DATAFLAGparams.matchObjPosSMAP_IMG==1
                        I_alias_sense_temp_7d = I_alias_sense_7d*0;
                        delta = DATAFLAGparams.matchObjPosDelta(ind_chunk);
                        if strcmpi(DWIparams.scan_mode,'3D')
                            [ny1,nx1,nz1,nchunk1,nc1,ndw1,nori1,nshot1] = size(I_alias_sense_7d);
                        else
                            [ny1,nx1,nz1,nc1,ndw1,nori1,nshot1] = size(I_alias_sense_7d);
                        end
                        %* nz1 and nchunk1 are always 1.
                        for ind1=1:nc1
                            for ind2 = 1:ndw1
                                for ind3 = 1:nori1
                                    for ind4 = 1:nshot1
                                        if strcmpi(DWIparams.fat_shift_dir,'P')
                                            if strcmpi(DWIparams.scan_mode,'3D')
                                                I_alias_sense_temp_7d(:,:,:,:,ind1,ind2,ind3,ind4) = ...
                                                    circshift(squeeze( ...
                                                    I_alias_sense_7d(:,:,:,:,ind1,ind2,ind3,ind4)), ...
                                                    [-round(delta),0,0,0]);
                                            else
                                                I_alias_sense_temp_7d(:,:,:,ind1,ind2,ind3,ind4) = ...
                                                    circshift(squeeze( ...
                                                    I_alias_sense_7d(:,:,:,ind1,ind2,ind3,ind4)), ...
                                                    [-round(delta),0,0]);
                                            end
                                        elseif strcmpi(DWIparams.fat_shift_dir,'L')
                                            if strcmpi(DWIparams.scan_mode,'3D')
                                                I_alias_sense_temp_7d(:,:,:,:,ind1,ind2,ind3,ind4) = ...
                                                    circshift(squeeze( ...
                                                    I_alias_sense_7d(:,:,:,:,ind1,ind2,ind3,ind4)), ...
                                                    [0,-round(delta),0,0]);
                                            else
                                                I_alias_sense_temp_7d(:,:,:,ind1,ind2,ind3,ind4) = ...
                                                    circshift(squeeze( ...
                                                    I_alias_sense_7d(:,:,:,ind1,ind2,ind3,ind4)), ...
                                                    [0,-round(delta),0]);
                                            end
                                        else
                                            error('s_recon_DWI2:main','Unknown DWIparams.fat_shift_dir')
                                        end
                                    end
                                end
                            end
                        end
                        I_alias_sense_7d = I_alias_sense_temp_7d;
                        clear  I_alias_sense_temp_7d
                        fprintf('  I_alias_sense_7d is shifted by -delta=[%d]\n',-round(delta))
                    end
                    
                    
                    % Get mask slice.
                    I_mask_img_m = I_mask_S_img_3d(:,:,ind_slice1);
                    
                    
                    %-------------------- START --------------------
                    % Erode mask to eliminate very high intensity noisy voxel
                    % in SMAP.
                    
                    % Take care of Volume Coil case (REFparams.filename is
                    % empty).
                    if ~isempty(REFparams.filename)
                        
                        % This is default.
                        %erodeval = SMAPparams.vox_erode+3;
                        erodeval = SMAPparams.vox_erode;
                        I_mask_img_temp_m = I_mask_img_m;
                        %I_mask_img_temp_m([1:erodeval,Ny-erodeval+1:Ny],:) = 0;
                        I_mask_img_temp_m(:,[1:erodeval,Nx-erodeval+1:Nx]) = 0;
                        I_mask_img_temp_m = imerode(I_mask_img_temp_m,strel('disk',erodeval-1));
                        I_mask_img_m = I_mask_img_temp_m;
                        
                        if strcmpi(GENparams.coilID,'SENSE-Breast-4')
                            %erodeval = SMAPparams.vox_erode+3;
                            erodeval = round(SMAPparams.vox_dilate/1.5);
                            I_mask_img_temp_m = I_mask_img_m;
                            I_mask_img_temp_m(:,[1,end]) = 0;
                            I_mask_img_temp_m([1,end],:) = 0;
                            I_mask_img_temp_m = imerode(I_mask_img_temp_m,strel('disk',erodeval));
                            I_mask_img_m = I_mask_img_temp_m;
                        end
                        clear  I_mask_img_temp_m  erodeval
                    end
                    %-------------------- END --------------------
                    
                    
                    % Generate aliased image data for the whole NSA (navg).
                    g_m = zeros(Ny/(R*Ns)*Nc*Ns*navg,nx);
                    for ind_avg = avg_v%1:navg
                        if isempty(flag_average_idx)
                            ind_avg1 = ind_avg; % use all navg
                        else
                            ind_avg1 = 1; % use single average
                        end
                        
                        % Generate aliased images for phase correction.
                        %rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1))+(Ny/(R*Ns)*Nc*Ns)*(ind_avg-1);
                        rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1))+(Ny/(R*Ns)*Nc*Ns)*(ind_avg1-1);
                        
                        % Change aliased image order for coil and shots to
                        % incorporate noise correlation.
                        if RECONFLAGparams.psiMethod==0 || RECONFLAGparams.psiMethod==1
                            %*** Method 0 or 1.
                            for ind_shot = 1:Ns
                                for ind_coil = 1:Nc
                                    ind = ind_coil+(ind_shot-1)*Nc;
                                    r1 = rstart_v(ind);
                                    r2 = r1+Ny/(R*Ns)-1;
                                    %g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                    %    :,ind_slice,ind_coil,ind_avg,ind_dw,ind_shot); % 2010.11.07
                                    if strcmpi(DWIparams.scan_mode,'3D')
                                        g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                            :,slice_idx,chunk_idx,ind_coil,ind_avg,ind_dw,ind_shot);
                                    else
                                        g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                            :,slice_idx,ind_coil,ind_avg,ind_dw,ind_shot);
                                    end
                                    
                                    if ind_coil==1
                                        rr1 = r1;
                                    end
                                    if ind_coil==Nc
                                        rr2 = r2;
                                    end
                                end
                                % Decorrelate g_m.
                                g_m(rr1:rr2,:) = PP * g_m(rr1:rr2,:);
                            end
                        elseif RECONFLAGparams.psiMethod==2
                            %*** Method 2.
                            for ind_coil = 1:Nc
                                for ind_shot = 1:Ns
                                    ind = ind_shot+(ind_coil-1)*Ns;
                                    r1 = rstart_v(ind);
                                    r2 = r1+Ny/(R*Ns)-1;
                                    %g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                    %    :,ind_slice,ind_coil,ind_avg,ind_dw,ind_shot); % 2010.11.07
                                    if strcmpi(DWIparams.scan_mode,'3D')
                                        g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                            :,slice_idx,chunk_idx,ind_coil,ind_avg,ind_dw,ind_shot);
                                    else
                                        g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                            :,slice_idx,ind_coil,ind_avg,ind_dw,ind_shot);
                                    end
                                end
                            end
                            % Decorrelate g_m.
                            rr1 = rstart_v(1);
                            rr2 = r2;
                            g_m(rr1:rr2,:) = PP * g_m(rr1:rr2,:);
                        end
                    end % for ind_avg
                    
                    
                    % Set max threshold.
                    if strcmpi(DWIparams.scan_mode,'3D')
                        val = squeeze(I_alias_sense_7d(:,:,slice_idx,chunk_idx,:,:,ind_dw,:));
                    else
                        val = squeeze(I_alias_sense_7d(:,:,slice_idx,:,:,ind_dw,:));
                    end
                    aliasedMax = max(abs(val(:)));  clear val
                    
                    % Output for each slice recon.
                    img_recon_m = zeros(Ny,nx);
                    %img_gfactor_m = zeros(Ny,nx);
                    %img_covar_m = zeros(Ny,nx);
                    
                    
                    % SENSE recon with phase correction. Loop through columns.
                    fprintf('    STARTING DWI recon ORI[%d], DW[%d], CHUNK[%d], SLICE[%d]\n', ...
                        dw_ori,ind_dw,ind_chunk,ind_slice)
                    
                    
                    % Get EPI factor.
                    EPI_FACTOR = Ny/(R*Ns);
                    
                    
                    %-----------------------------------
                    % Precalculation of encoding matrix
                    %-----------------------------------
                    
                    % Precalculation of E_m.
                    if RECONFLAGparams.precalc_encoding_matrix_flag==true
                        tic
                        E_c = cell(1,nx);
                        
                        for col = 1:nx
                            %if mod(col,50)==0
                            %    fprintf('. ')
                            %end
                            aliased_v = g_m(:,col);
                            
                            % Skip this column if there's little signal:
                            if (max(abs(aliased_v)) < 0.01*aliasedMax)
                                continue
                            end
                            if ~any(I_mask_img_m(:,col)==1)
                                continue
                            end
                            
                            % Pre-calculate encoding matrix.
                            % Method 1 and 2 work ok. Method 0 doesn' use PSI.
                            if RECONFLAGparams.psiMethod==0 || RECONFLAGparams.psiMethod==1
                                
                                %*** Method 0 or 1.
                                E_m = zeros(Ny/(R*Ns)*Nc*Ns*navg,Ny);
                                for ind_avg = avg_v%1:navg
                                    if isempty(flag_average_idx)
                                        ind_avg1 = ind_avg; % use all navg
                                    else
                                        ind_avg1 = 1; % use single average
                                    end
                                    %rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1))+(Ny/(R*Ns)*Nc*Ns)*(ind_avg-1);
                                    rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1))+(Ny/(R*Ns)*Nc*Ns)*(ind_avg1-1);
                                    for ind_shot = 1:Ns
                                        %eval(sprintf('P_m = P%d_c{ind_avg};',ind_shot))
                                        %eval(sprintf('P_v = P%d_c{ind_avg};',ind_shot))
                                        if strcmpi(DWIparams.scan_mode,'3D')
                                            phi_v = exp(1i*angle(squeeze(...
                                                nav_recon_6d(:,col,ind_slice,ind_chunk, ...
                                                ind_avg,ind_dw,ind_shot))));
                                        else
                                            phi_v = exp(1i*angle(squeeze(...
                                                nav_recon_6d(:,col,ind_slice, ...
                                                ind_avg,ind_dw,ind_shot))));
                                        end
                                        if (RECONFLAGparams.phaseCorr==1) && (dw_ori~=0)
                                            P_v = conj(phi_v);
                                        else
                                            P_v = ones(length(phi_v),1);
                                        end
                                        
                                        A_m = double(A_c{ind_shot});
                                        for ind_coil = 1:Nc
                                            ind = ind_coil+(ind_shot-1)*Nc;
                                            r1 = rstart_v(ind);
                                            r2 = r1+Ny/(R*Ns)-1;
                                            
                                            % Original
                                            %E_m(r1:r2,:) = A_m*diag(I_smap_img_4d(:,col,ind_slice,ind_coil))*P_m;
                                            
                                            % Modified
                                            %s_p_v = double(I_smap_img_4d(:,col,ind_slice1,ind_coil).* P_v);
                                            s_p_v = double(I_smap_img_4d(:,col,ind_slice2,ind_coil).* P_v);
                                            ASP_m = A_m;
                                            for indr = 1:EPI_FACTOR
                                                col_begin = indr;
                                                a_v = A_m(indr,:);
                                                ASP_m(indr,col_begin:EPI_FACTOR:end) = ...
                                                    a_v(col_begin:EPI_FACTOR:end).*...
                                                    s_p_v(col_begin:EPI_FACTOR:end).';
                                            end
                                            E_m(r1:r2,:) = ASP_m;
                                            
                                            % Reserve indices for PP.
                                            if ind_coil==1
                                                r00 = r1;
                                            end
                                            if ind_coil==Nc
                                                r11 = r2;
                                            end
                                        end
                                        
                                        % Decorrelate E_m.
                                        E_m(r00:r11,:) = sparse(PP)*sparse(E_m(r00:r11,:));
                                    end
                                end % for ind_avg
                                
                            elseif RECONFLAGparams.psiMethod==2
                                
                                %*** Method 2.
                                E_m = zeros(Ny/(R*Ns)*Nc*Ns*navg,Ny);
                                for ind_avg = avg_v%1:navg
                                    if isempty(flag_average_idx)
                                        ind_avg1 = ind_avg; % use all navg
                                    else
                                        ind_avg1 = 1; % use single average
                                    end
                                    %rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1))+(Ny/(R*Ns)*Nc*Ns)*(ind_avg-1);
                                    rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1))+(Ny/(R*Ns)*Nc*Ns)*(ind_avg1-1);
                                    for ind_coil = 1:Nc
                                        for ind_shot = 1:Ns
                                            %eval(sprintf('P_m = P%d_c{ind_avg};',ind_shot))
                                            if strcmpi(DWIparams.scan_mode,'3D')
                                                phi_v = exp(1i*angle(squeeze(...
                                                    nav_recon_6d(:,col,ind_slice,ind_chunk, ...
                                                    ind_avg,ind_dw,ind_shot))));
                                            else
                                                phi_v = exp(1i*angle(squeeze(...
                                                    nav_recon_6d(:,col,ind_slice, ...
                                                    ind_avg,ind_dw,ind_shot))));
                                            end
                                            if (RECONFLAGparams.phaseCorr==1) && (dw_ori~=0)
                                                P_v = conj(phi_v);
                                            else
                                                P_v = ones(length(phi_v),1);
                                            end
                                            
                                            A_m = A_c{ind_shot};
                                            ind = ind_shot + (ind_coil-1)*Ns;
                                            r1 = rstart_v(ind);
                                            r2 = r1+Ny/(R*Ns)-1;
                                            
                                            % Original.
                                            %E_m(r1:r2,:) = A_m*diag(I_smap_img_4d(:,col,ind_slice,ind_coil))*P_m;
                                            
                                            % Modified.
                                            %s_p_v = double(I_smap_img_4d(:,col,ind_slice1,ind_coil).*P_v);
                                            s_p_v = double(I_smap_img_4d(:,col,ind_slice2,ind_coil).*P_v);
                                            ASP_m = A_m;
                                            for indr = 1:EPI_FACTOR
                                                col_begin = indr;
                                                a_v = A_m(indr,:);
                                                ASP_m(indr,col_begin:EPI_FACTOR:end) = ...
                                                    a_v(col_begin:EPI_FACTOR:end).*...
                                                    s_p_v(col_begin:EPI_FACTOR:end).';
                                            end
                                            E_m(r1:r2,:) = ASP_m;
                                        end
                                    end
                                    
                                    % Decorrelate E_m.
                                    %* PP is multiplied repeatedly over averages
                                    E_m(rstart_v(1):r2,:) = sparse(PP)*sparse(E_m(rstart_v(1):r2,:));
                                    
                                end % for ind_avg
                            end % if RECONFLAGparams.psiMethod
                            
                            % Reserve encoding matrix.
                            E_c{col} = E_m;
                            
                        end % for col
                        t0 = toc;
                        fprintf('\n')
                        fprintf('    %f sec took for reserving all E\n',t0)
                    end % if RECONFLAGparams.precalc_encoding_matrix_flag
                    
                    
                    
                    %--------------------------------------------------------
                    % Inversion of E, precalc or columnwise w/ or w/o sub-E
                    %--------------------------------------------------------
                    
                    % Determine sub-encoding matrix size. Force the column size
                    % of encoding matrix to be Ny when precalculation is used
                    % or sub-encoding matrix is not used.
                    % Precalculation with sub-encoding matrix (e.g., R*Ns)
                    % didn't make much difference.
                    if RECONFLAGparams.precalc_encoding_matrix_flag == true
                        matrix_siz = Ny;
                    else
                        if RECONFLAGparams.sub_encoding_matrix_flag==true
                            sub_enc_idx = ...
                                RECONFLAGparams.sub_encoding_matrix_siz_idx;
                            switch sub_enc_idx
                                case 0
                                    matrix_siz = Ny;
                                case 1
                                    matrix_siz = R;
                                case 2
                                    matrix_siz = Ns;
                                case 3
                                    matrix_siz = R*Ns;
                                otherwise
                                    error('s_recon_DWI__precalc:main',...
                                        'Unknown sub_encoding_matrix_siz_idx')
                            end
                        else
                            matrix_siz = Ny;
                        end
                    end
                    
                    % Invert encoding matrix and recon each image column.
                    tic
                    for col = 1:nx
                        aliased_v = g_m(:,col);
                        
                        % Skip this column if there's little signal:
                        if (max(abs(aliased_v)) < 0.01*aliasedMax)
                            continue
                        end
                        if ~any(I_mask_img_m(:,col)==1)
                            continue
                        end
                        
                        % Retrieve or calculate encoding matrix.
                        if RECONFLAGparams.precalc_encoding_matrix_flag == true
                            E_m = E_c{col};
                        else
                            % Calculate encoding matrix, if pre-calculation is
                            % not used.
                            % Method 1 and 2 work ok. Method 0 doesn' use PSI.
                            if RECONFLAGparams.psiMethod==0 || ...
                                    RECONFLAGparams.psiMethod==1
                                
                                % Method 0 or 1.
                                E_m = zeros(Ny/(R*Ns)*Nc*Ns*navg,Ny);
                                for ind_avg = avg_v%1:navg
                                    if isempty(flag_average_idx)
                                        ind_avg1 = ind_avg; % use all navg
                                    else
                                        ind_avg1 = 1; % use single average
                                    end
                                    %rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1))+(Ny/(R*Ns)*Nc*Ns)*(ind_avg-1);
                                    rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1))+(Ny/(R*Ns)*Nc*Ns)*(ind_avg1-1);
                                    for ind_shot = 1:Ns
                                        %eval(sprintf('P_m = P%d_c{ind_avg};',ind_shot))
                                        %eval(sprintf('P_v = P%d_c{ind_avg};',ind_shot))
                                        if strcmpi(DWIparams.scan_mode,'3D')
                                            phi_v = exp(1i*angle(squeeze(...
                                                nav_recon_6d(:,col,ind_slice,ind_chunk, ...
                                                ind_avg,ind_dw,ind_shot))));
                                            
                                            % Calculate relative phase
                                            % error to ref slice (kz).
%                                             if ind_slice==ref
%                                                 phi_m = exp(1i*angle(squeeze(...
%                                                     nav_recon_6d(:,:,ind_slice,ind_chunk, ...
%                                                     ind_avg,ind_dw,ind_shot))));
%                                                 eval(sprintf('phi%d_m = phi_m;',ind_shot))
%                                             else
%                                                 eval(sprintf('phi_m = phi%d_m;',ind_shot))
%                                                 phi0_v = phi_m(:,col);
%                                                 phi_v = phi0_v .*conj(phi_v); % phase difference to phi0_v at reference kz
%                                             end
                                        else
                                            phi_v = exp(1i*angle(squeeze(...
                                                nav_recon_6d(:,col,ind_slice, ...
                                                ind_avg,ind_dw,ind_shot))));
                                        end
                                        if (RECONFLAGparams.phaseCorr==1) && (dw_ori~=0)
                                            P_v = conj(phi_v);
                                        else
                                            P_v = ones(length(phi_v),1);
                                        end
                                        
                                        A_m = A_c{ind_shot};
                                        for ind_coil = 1:Nc
                                            ind = ind_coil+(ind_shot-1)*Nc;
                                            r1 = rstart_v(ind);
                                            r2 = r1+Ny/(R*Ns)-1;
                                            
                                            % Original
                                            %E_m(r1:r2,:) = A_m*diag(I_smap_img_4d(:,col,ind_slice,ind_coil))*P_m;
                                            
                                            % Modified
                                            %s_p_v = double(I_smap_img_4d(:,col,ind_slice1,ind_coil).*P_v);
                                            s_p_v = double(I_smap_img_4d(:,col,ind_slice2,ind_coil).*P_v);
                                            ASP_m = A_m;
                                            for indr = 1:EPI_FACTOR
                                                col_begin = indr;
                                                a_v = A_m(indr,:);
                                                ASP_m(indr,col_begin:EPI_FACTOR:end) = ...
                                                    a_v(col_begin:EPI_FACTOR:end).*...
                                                    s_p_v(col_begin:EPI_FACTOR:end).';
                                            end
                                            E_m(r1:r2,:) = ASP_m;
                                            
                                            % Reserve indices for PP.
                                            if ind_coil==1
                                                r00 = r1;
                                            end
                                            if ind_coil==Nc
                                                r11 = r2;
                                            end
                                        end
                                        
                                        % Decorrelate E_m.
                                        E_m(r00:r11,:) = sparse(PP)*sparse(E_m(r00:r11,:));
                                    end
                                end % for ind_avg
                                
                            elseif RECONFLAGparams.psiMethod==2
                                
                                % Method 2.
                                E_m = zeros(Ny/(R*Ns)*Nc*Ns*navg,Ny);
                                for ind_avg = avg_v%1:navg
                                    if isempty(flag_average_idx)
                                        ind_avg1 = ind_avg; % use all navg
                                    else
                                        ind_avg1 = 1; % use single average
                                    end
                                    %rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1))+(Ny/(R*Ns)*Nc*Ns)*(ind_avg-1);
                                    rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1))+(Ny/(R*Ns)*Nc*Ns)*(ind_avg1-1);
                                    for ind_coil = 1:Nc
                                        for ind_shot = 1:Ns
                                            %eval(sprintf('P_m = P%d_c{ind_avg};',ind_shot))
                                            if strcmpi(DWIparams.scan_mode,'3D')
                                                phi_v = exp(1i*angle(squeeze(...
                                                    nav_recon_6d(:,col,ind_slice,ind_chunk, ...
                                                    ind_avg,ind_dw,ind_shot))));
                                            else
                                                phi_v = exp(1i*angle(squeeze(...
                                                    nav_recon_6d(:,col,ind_slice, ...
                                                    ind_avg,ind_dw,ind_shot))));
                                            end
                                            if (RECONFLAGparams.phaseCorr==1) && (dw_ori~=0)
                                                P_v = conj(phi_v);
                                            else
                                                P_v = ones(length(phi_v),1);
                                            end
                                            
                                            A_m = A_c{ind_shot};
                                            ind = ind_shot + (ind_coil-1)*Ns;
                                            r1 = rstart_v(ind);
                                            r2 = r1+Ny/(R*Ns)-1;
                                            
                                            % Original.
                                            %E_m(r1:r2,:) = A_m*diag(I_smap_img_4d(:,col,ind_slice,ind_coil))*P_m;
                                            
                                            % Modified.
                                            %s_p_v = double(I_smap_img_4d(:,col,ind_slice1,ind_coil).*P_v);
                                            s_p_v = double(I_smap_img_4d(:,col,ind_slice2,ind_coil).*P_v);
                                            ASP_m = A_m;
                                            for indr = 1:EPI_FACTOR
                                                col_begin = indr;
                                                a_v = A_m(indr,:);
                                                ASP_m(indr,col_begin:EPI_FACTOR:end) = ...
                                                    a_v(col_begin:EPI_FACTOR:end).*...
                                                    s_p_v(col_begin:EPI_FACTOR:end).';
                                            end
                                            E_m(r1:r2,:) = ASP_m;
                                        end
                                    end
                                    
                                    % Decorrelate E_m.
                                    %* PP is multiplied repeatedly over averages
                                    E_m(rstart_v(1):r2,:) = sparse(PP)*sparse(E_m(rstart_v(1):r2,:));
                                    
                                end % for ind_avg
                            end % if RECONFLAGparams.psiMethod
                        end % if RECONFLAGparams.precalc_encoding_matrix_flag
                        
                        
                        % Invert sub- or full-encoding matrix.
                        if matrix_siz < Ny
                            % Calculate unfolded and phase-corrected image for
                            % part of each column.
                            leap_siz = size(E_m,2)/matrix_siz;
                            for ind_row = 1:leap_siz
                                col_pos_v = ind_row:leap_siz:size(E_m,2);
                                E_sub_m = E_m(:,col_pos_v);
                                cImage_sub_v = E_sub_m \ aliased_v;
                                cImage_sub_v(isnan(cImage_sub_v)) = 0;
                                
                                % Save recon image
                                img_recon_m(col_pos_v,col) = cImage_sub_v;
                                %img_gfactor_m(col_pos_v,col) = diag(sqrt( ...
                                %    pinv(full(E_sub_m'*E_sub_m)) .* (full(E_sub_m'*E_sub_m)) ));
                            end
                        else
                            % Calculate unfolded and phase-corrected image for
                            % full column.
                            A = E_m'*E_m;
                            b = E_m'*aliased_v;
                            S = sparse(A);
                            cImage_v = S \ b;
                            cImage_v(isnan(cImage_v))=0;
                            
                            % Save recon image
                            img_recon_m(:,col) = cImage_v;
                            %img_gfactor_m(:,col) = diag(sqrt( ...
                            %    pinv(full(E_m'*E_m)) .* (full(E_m'*E_m)) ));
                        end
                    end % for col
                    t1 = toc;
                    fprintf('    Finished unfolding DWI in [%.3f]sec\n',t1)
                    
                    
                    %----- TEMP -----
                    %sendEmail(sprintf('slice[%d],chunk[%d],ori[%d] recon in [%.3f]sec at [%s]', ...
                    %    ind_slice,ind_chunk,dw_ori,t1,DWIparams.filename),' ')
                    %----- TEMP -----
                    
                    
                    % Mask recon image.
                    m = img_recon_m.*I_mask_img_m;
                    %g = img_gfactor_m.*I_mask_img_m;
                    %showimage(gray,abs(m.*imerode(I_mask_img_m,strel('disk',14))), ...
                    %    angle(m.*imerode(I_mask_img_m,strel('disk',14))),sprintf('slice %d',ind_slice))
                    
                    clear  img_recon_m  img_gfactor_m  E_m  cImage_v  A  b  S
                    
                    % Keep output.
                    if strcmpi(DWIparams.scan_mode,'3D')
                        img_recon_4d(:,:,ind_slice,ind_chunk,ind_dw) = m;
                        %img_gfactor_4d(:,:,ind_slice,ind_chunk,ind_dw) = g;
                    else
                        img_recon_4d(:,:,ind_slice,ind_dw) = m;
                        %img_gfactor_4d(:,:,ind_slice,ind_dw) = g;
                    end
                    %clear  m  g
                    
                    if ind_slice==4
                        eval(sprintf('save  temp_slice%.2d_chunk%.2d_ori%.2d  m',...
                            ind_slice,ind_chunk,dw_ori))
                    end
                    
                    % Test show image.
                    %showimage(jet,abs(m),angle(m),sprintf('SLICE[%d],DW\\_ORI[%d]',ind_slice,dw_ori))
                    
                    % TEMP.
%                     m_covar = img_covar_m.*I_mask_img_m;
%                     if RECONFLAGparams.psiMethod==0
%                         m0=m;
%                         m_covar0=m_covar;
%                         snr0 = abs(m0)./sqrt(abs(m_covar0));
%                         showimage(jet,abs(m0),abs(m_covar0),'img0','covar0')
%                         showimage(jet,abs(m0),snr0,'img0','snr0')
%                     else
%                         m1 = m;
%                         m_covar1 = m_covar;
%                         snr1 = abs(m1)./sqrt(abs(m_covar1));
%                         showimage(jet,abs(m1),abs(m_covar1),'img1','covar1')
%                         showimage(jet,abs(m1),snr1,'img1','snr1')
%                     end
%                     
%                     %
%                     [mean(abs(m0(bw(:)))),mean(abs(m1(bw(:))))]
%                     [mean(abs(m_covar0(bw(:)))),mean(abs(m_covar1(bw(:))))]
%                     sqrt([mean(abs(m_covar0(bw(:)))),mean(abs(m_covar1(bw(:))))])
%                     [mean(abs(snr0(bw(:)))),mean(abs(snr1(bw(:))))]
%                     [std(real(noi_dwi_I_v)),std(imag(noi_dwi_I_v)),mean(abs(noi_dwi_I_v))/sqrt(pi/2)]
                    
                end % for ind_slice
                
                % IFT along kz for 3D acquisition.
                img_recon_temp_4d = img_recon_4d; %------TEMP------
                if strcmpi(DWIparams.scan_mode,'3D')
                    for ind1 = 1:ndw
%                         img_recon_4d(:,:,:,ind_chunk,ind1) = ...
%                             fftshift(ift1(squeeze(img_recon_4d(:,:,:,ind_chunk,ind1)),3),3);
                        img_recon_4d(:,:,:,ind_chunk,ind1) = ...
                            ift1(squeeze(img_recon_4d(:,:,:,ind_chunk,ind1)),3);
                    end
                end
                fprintf('\n*** IFT along kz is performed for chunk[%d],ori[%d] ***\n', ...
                    ind_chunk,dw_ori)
                
            end % for ind_chunk
        end % for ind_dw
        clear  I_alias_sense_7d
        
        
        % Save recon data.
        cd(saveReconDir_s)
        if isempty(flag_average_idx)
            fname_s = sprintf('%s_avg%.2d_ori%.2d__R%dS%d', ...
                fname_img_head_s,navg,dw_ori,R,Ns);
        else
            fname_s = sprintf('%s_avg%.2d_ori%.2d__R%dS%d_avgIdx%d', ...
                fname_img_head_s,navg,dw_ori,R,Ns,flag_average_idx);
        end
        %eval(sprintf('save  %s  img_recon_4d  img_gfactor_4d',fname_s))
        eval(sprintf('save  %s  img_recon_4d',fname_s))
        
        fprintf('\n')
        fprintf('    [%s] saved at \n',fname_s)
        fprintf('    [%s]\n\n',saveReconDir_s)
        clear  img_recon*  mask  fname_s  m  g  I_mask_img_m  P_* phi*
        
    end % for dw_ori
    
    % Clear unused data.
    clear  I_smap_*  I_mask_*  E_*
    
end % if (RECONFLAGparams.phaseCorr)==1 && (BATCHparams.recon_data==1)

% Pack memory.
pack

fprintf('\n\n')









