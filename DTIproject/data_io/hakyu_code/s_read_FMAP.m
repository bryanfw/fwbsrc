
%[s_read_FMAP] reads, resize and save fieldmap in PAR/REC format.
%
%
% Last modified
% 2010.02.11.
% 2010.02.17. saveReconDir_s -> saveDataDir_s to separate DW data and recon data.
% 2010.08.09
%   This is generated from [FMAP.m].
% 2010.08.18.
%   Pack memory.
% 2010.08.20.
%   Modify for matching object size of FMAP (resizing) and IMG.
% 2010.08.22.
%   Check the resize step of FMAP, Eq.[1] and [2].
% 2010.09.11.
%   Modify R = RECONparams.ENCima.sense_factors(2); to R =
%   DWIparams.SENSE_FACTOR; to eliminate non-integral sense factor.
% 2010.09.20.
%   Adjust slice_orientation.
% 2010.12.09.
%   Adjust slice_orientation of FMAP data for resizing and recon for DWI
%   data.
% 2012.05.16.
%   Don't read FMAP when FMAPparams.filename is empty.
%
% Ha-Kyu



%% Read, resize and save PAR/REC fieldmap

fprintf('Read, resize and save PAR/REC fieldmap.\n')

% Read data.
if DATAFLAGparams.genFMAP == 1
    
    if ~isempty(FMAPparams.filename)
        cd(loadParDir_s)
        [data,info] = f_read_parrec(FMAPparams.filename); % [240x240] and recon into [1.0x1.0]mm
        clear  info
        fprintf('    [%s] is read as fieldmap.\n',FMAPparams.filename)
    else
        fprintf('    [%s] is read as fieldmap.\n',FMAPparams.filename)
        break
    end
    
    % Get fieldmap.
    if length(size(data))~=8
        fprintf('\n')
        fprintf('    Size of FMAP data is not what is expected. Then skip reading FMAP.\n')        
        FMAPparams.filename = [];
        fprintf('    FMAPparams.filename is set to []\n')
        fprintf('\n\n')
        return
    end
    magn_3d = squeeze(data(:,:,:,1,1,1, 1,1));
    fmap_3d = squeeze(data(:,:,:,1,1,1, 2,2));
    [ny,nx,nsl] = size(magn_3d);
    clear  data

    
    % Get FMAP image size parameter.
    % f_N_res_p(m), f_zeropad_vox_p(m), f_cut_vox_p(m) are all along the
    % P_ORI(M_ORI) of the main DWI data. Then FMAP data must permute,
    % fliplr or flipud based on 'DWIparams.fold_over_dir'.
    
    % Image object size matching along P_ORI can be done as,
    %   ovs_fac(2)*FOV(1) = ovs_fac(1)*FOV(2),
    % here, (1) and (2) are FMAP and IMG, respectively. Then,
    %   FOV(1) = ovs_fac(1)*FOV(2)/ovs_fac(2),
    % and ovs_fac(1)=1, then
    %   FOV(1) = N(1)*deltaY(1) = N(2)*deltaY(2)/ovs_fac(2),
    % and only matrix index matters in image matrix, then
    %   N(1) = N(2) / ovs_fac(2),
    % and because of epi factor, ns and R,
    %   N(1) = E*ns*R / ovs_fac(2) --- Eq.[1]
    %
    % Similarly, along M_ORI, can be,
    %   N(1) = N(2) / ovs_fac(2) --- Eq.[2]
    
    %if strcmpi(DWIparams.fold_over_dir,'AP')   
    
        %* P_ORI - always matrix row direction in DWI.
        E = DWIparams.EPI_FACTOR(1);
        R = DWIparams.SENSE_FACTOR;
        ns = RECONparams.ACQREC.number_of_shots;
        ovs_fac = RECONparams.ENCima.oversample_factors(2);
        
        f_N_res_p = (E*ns*R)/ovs_fac; % Eq.[1]
        
        if mod((E*ns*R) - ceil(f_N_res_p),2)==0
            f_N_res_p = ceil(f_N_res_p);
        else
            f_N_res_p = floor(f_N_res_p);
        end
        f_zeropad_vox_p = [];
        f_cut_vox_p = [];
        if f_N_res_p < (E*ns*R)
            f_zeropad_vox_p = ((E*ns*R) - f_N_res_p)/2;
        else
            f_cut_vox_p = (f_N_res_p - (E*ns*R))/2;
        end
        
        %* M_ORI - always matrix column direction in DWI.
        N = RECONparams.ENCima.oversample_resolutions(1);
        ovs_fac = RECONparams.ENCima.oversample_factors(1);
        
        f_N_res_m = N/ovs_fac; % Eq.[2]
        
        if mod(N - ceil(f_N_res_m),2)==0
            f_N_res_m = ceil(f_N_res_m);
        else
            f_N_res_m = floor(f_N_res_m);
        end
        f_zeropad_vox_m = [];
        f_cut_vox_m = [];
        if f_N_res_m < N
            f_zeropad_vox_m = (N - f_N_res_m)/2;
        else
            f_cut_vox_m = (f_N_res_m - N)/2;
        end        
                
%     elseif strcmpi(DWIparams.fold_over_dir,'RL')
%         
%         %* P_ORI.
%         N = RECONparams.ENCima.oversample_resolutions(1);
%         ovs_fac = RECONparams.ENCima.oversample_factors(1);
%         
%         f_N_res_p = N/ovs_fac;
%         
%         if mod(N - ceil(f_N_res_p),2)==0
%             f_N_res_p = ceil(f_N_res_p);
%         else
%             f_N_res_p = floor(f_N_res_p);
%         end
%         f_zeropad_vox_p = [];
%         f_cut_vox_p = [];
%         if f_N_res_p < N
%             f_zeropad_vox_p = (N - f_N_res_p)/2;
%         else
%             f_cut_vox_p = (f_N_res_p - N)/2;
%         end
%         
%                         
%         %* M_ORI.
%         E = DWIparams.EPI_FACTOR(1);
%         R = DWIparams.SENSE_FACTOR;
%         ns = RECONparams.ACQREC.number_of_shots;
%         ovs_fac = RECONparams.ENCima.oversample_factors(2);
%         
%         f_N_res_m = (E*ns*R)/ovs_fac;
%         
%         if mod(E*ns*R - ceil(f_N_res_m),2)==0
%             f_N_res_m = ceil(f_N_res_m);
%         else
%             f_N_res_m = floor(f_N_res_m);
%         end
%         f_zeropad_vox_m = [];
%         f_cut_vox_m = [];
%         if f_N_res_m < (E*ns*R)
%             f_zeropad_vox_m = ((E*ns*R) - f_N_res_m)/2;
%         else
%             f_cut_vox_m = (f_N_res_m - (E*ns*R))/2;
%         end
%         
%         %* Adjust FMAP object position for reconstruction variable with
%         %* the main DWI data foldover and fatshift.
%         if strcmpi(DWIparams.fat_shift_dir,'L')
%             for ind_slice = 1:size(magn_3d,3)
%                 magn_3d(:,:,ind_slice) = permute(flipud(magn_3d(:,:,ind_slice)),[2,1]);
%                 fmap_3d(:,:,ind_slice) = permute(flipud(fmap_3d(:,:,ind_slice)),[2,1]);
%             end
%         elseif strcmpi(DWIparams.fat_shift_dir,'R')
%             for ind_slice = 1:size(magn_3d,3)
%                 magn_3d(:,:,ind_slice) = permute(fliplr(magn_3d(:,:,ind_slice)),[2,1]);
%                 fmap_3d(:,:,ind_slice) = permute(fliplr(fmap_3d(:,:,ind_slice)),[2,1]);
%             end
%         end            
%     end


    % Adjust FMAP object position for resizing and recon with the main DWI
    % data foldover and fatshift direction.
    if (strcmpi(DWIparams.slice_orientation,'transverse') && ...
            strcmpi(DWIparams.fat_shift_dir,'L')) || ...
            (strcmpi(DWIparams.slice_orientation,'coronal') && ...
            strcmpi(DWIparams.fat_shift_dir,'L')) || ...
            (strcmpi(DWIparams.slice_orientation,'sagittal') && ...
            strcmpi(DWIparams.fat_shift_dir,'P'))
            
        for ind_slice = 1:size(magn_3d,3)
            magn_3d(:,:,ind_slice) = permute(flipud(magn_3d(:,:,ind_slice)),[2,1]);
            fmap_3d(:,:,ind_slice) = permute(flipud(fmap_3d(:,:,ind_slice)),[2,1]);
        end
        
    elseif (strcmpi(DWIparams.slice_orientation,'transverse') && ...
            strcmpi(DWIparams.fat_shift_dir,'R')) || ...
            (strcmpi(DWIparams.slice_orientation,'coronal') && ...
            strcmpi(DWIparams.fat_shift_dir,'R')) || ...
            (strcmpi(DWIparams.slice_orientation,'sagittal') && ...
            strcmpi(DWIparams.fat_shift_dir,'A'))
        
        for ind_slice = 1:size(magn_3d,3)
            magn_3d(:,:,ind_slice) = permute(fliplr(magn_3d(:,:,ind_slice)),[2,1]);
            fmap_3d(:,:,ind_slice) = permute(fliplr(fmap_3d(:,:,ind_slice)),[2,1]);
        end
        
    elseif strcmpi(DWIparams.slice_orientation,'transverse') && ...
            strcmpi(DWIparams.fat_shift_dir,'P')
        
        % do nothing
        
    elseif strcmpi(DWIparams.slice_orientation,'transverse') && ...
            strcmpi(DWIparams.fat_shift_dir,'A')
        
        for ind_slice = 1:size(magn_3d,3)
            magn_3d(:,:,ind_slice) = fliplr(flipud(magn_3d(:,:,ind_slice)));
            fmap_3d(:,:,ind_slice) = fliplr(flipud(fmap_3d(:,:,ind_slice)));
        end
        
    else
        error('s_read_FMAP:main','Unknown DWI slice_orientation and fat_shift_dir')
    end
    [ny,nx,nsl] = size(magn_3d);
    
        
    % Resize fieldmap. This is final resize.
    magn_res_3d = zeros(RECONparams.s_N_final_p, RECONparams.s_N_final_m, nsl);
    fmap_res_3d = zeros(RECONparams.s_N_final_p, RECONparams.s_N_final_m, nsl);
    
    for ind_slice = 1:nsl
        m1 = magn_3d(:,:,ind_slice);
        m2 = fmap_3d(:,:,ind_slice);
        
        % Resize along P_ORI.
        m11 = imresize(m1,[f_N_res_p,f_N_res_m]);
        m21 = imresize(m2,[f_N_res_p,f_N_res_m]);
        
        % Zeropad or cut.
        if ~isempty(f_zeropad_vox_p)
            m11 = padarray(m11,[f_zeropad_vox_p,0],0,'both');
            m21 = padarray(m21,[f_zeropad_vox_p,0],0,'both');
        end
        if ~isempty(f_cut_vox_p)
            m11 = m11(f_cut_vox_p+1:end-f_cut_vox_p,:);
            m21 = m21(f_cut_vox_p+1:end-f_cut_vox_p,:);
        end
        if ~isempty(f_zeropad_vox_m)
            m11 = padarray(m11,[0,f_zeropad_vox_m],0,'both');
            m21 = padarray(m21,[0,f_zeropad_vox_m],0,'both');
        end
        if ~isempty(f_cut_vox_m)
            m11 = m11(:,f_cut_vox_m+1:end-f_cut_vox_m);
            m21 = m21(:,f_cut_vox_m+1:end-f_cut_vox_m);
        end
        
        % Keep resized fmap.
        magn_res_3d(:,:,ind_slice) = m11;
        fmap_res_3d(:,:,ind_slice) = m21;
    end
    fprintf('    fieldmap resized.\n')
    
    
    % Generate mask.
    mask_3d = zeros(size(magn_res_3d));
    for ind_slice = 1:nsl
        m = magn_res_3d(:,:,ind_slice);
        bw = m/max(m(:)) > graythresh(m)*0.1;
        bw = imdilate(bw,strel('disk',FMAPparams.vox_dilate),0);
        mask_3d(:,:,ind_slice) = bw;
    end
    fprintf('    mask generated.\n')
    
    
    % Output.
    I_fmap = struct('fmap',cat(4,magn_res_3d,fmap_res_3d),'mask',mask_3d);
    fprintf('    Resized fieldmap, F, generated.\n')
    
    
    % Save fieldmap.
    cd(saveDataDir_s)
    save  I_fmap  I_fmap
    fprintf('    fieldmap is saved.\n')
    
    
    % Clear data.
    clear  magn_3d fmap_3d  magn_res_3d  fmap_res_3d  mask_3d  bw m
    clear  m1  m2  m11  m21  X*  Y*  ovs_fac_m  Ny_ref_new  acq_resol
end

% Pack memory.
pack

fprintf('\n\n')



%% END








