
%[s_reslice_map_3d_v3] reslices SMAP and FMAP to match slice positions of DWI
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
% 2012.07.05, 2012.07.06
%   This version [s_reslice_map_v2.m] uses FOV and Stack Offc. information
%   in reslicing the SMAP.
%
%   Remove if statement checking the existence of [I_smap_img_4d.mat] file,
%   because DWIparams.multichunk related parameters must be saved.
% 2012.07.06.
%   [s_reslice_map_v2.m] is renamed to [s_reslice_map_3d.m].
%   [s_reslice_map_3d_v2.m] is generated.
% 2012.07.18.
%   Save I_smap_img_4d for each chunk in this version.
%
% Ha-Kyu



%% Resample slices of SMAP and FMAP for NAV and IMG recon

fprintf('Resample slices of SMAP and FMAP for NAV and IMG recon\n')



%% Check if SMAP is already resliced for IMG and NAV

cd(saveDataDir_s)
%if ~(exist('I_smap_img_4d.mat','file')==2) && ...
%        ~(exist('I_smap_nav_4d.mat','file')==2) && ...
%        ~isempty(REFparams.filename)
% if ~(exist('I_smap_img_4d.mat','file')==2) && ...
%         ~(exist('I_smap_img_sl01_4d.mat','file')) && ...
%         ~(isempty(REFparams.filename))
if ~(exist('I_smap_img_4d.mat','file')==2) && ...
        ~(exist('I_smap_img_chunk01_4d.mat','file')) && ...
        ~(isempty(REFparams.filename))
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
    
    % Keep nSLICE as shown in ExamCard. This will set ovs_fac to 1.
    nSLICE = DWIparams.nSLICE;
    
    
    
    % ---------- New routine for reslicing SMAP ----------
    
    % REF slice FOV and Offc.
    if strcmpi(REFparams.slice_orientation,'transverse')
        ref_fov = REFparams.FOV.FH;
        ref_offc = REFparams.Offc.FH;
    elseif strcmpi(REFparams.slice_orientation,'coronal')
        ref_fov = REFparams.FOV.AP;
        ref_offc = REFparams.Offc.AP;
    else
        ref_fov = REFparams.FOV.RL;
        ref_offc = REFparams.Offc.RL;
    end
    ref_th = REFparams.slice_thickness;
    ref_ovs_fac = snz / REFparams.slices; % ovs_fac along kz due to 3D acq.
    ref_fov_ovs = ref_fov * ref_ovs_fac;
    
    %* REF slice center.
    %** Slice 1 to end is from +ve to -ve coord.
    ref_mm_v = linspace(ref_fov_ovs/2 - ref_th/2,-ref_fov_ovs/2 + ref_th/2,snz);
    ref_mm_v = ref_mm_v + ref_offc; % REF slice position in gradient coord
        
    %* Check slice thickness.
    val = unique(abs(diff(ref_mm_v)));
    if val~=ref_th
        error('s_reslice_map_v2:main', ...
            'SMAP slice thickness must be the same as slice center locations')
    end
    
    
    % DWI slice FOV and Offc.
    if strcmpi(DWIparams.slice_orientation,'transverse')
        dwi_fov = DWIparams.FOV.FH;
        dwi_offc = DWIparams.Offc.FH;
    elseif strcmpi(DWIparams.slice_orientation,'coronal')
        dwi_fov = DWIparams.FOV.AP;
        dwi_offc = DWIparams.Offc.AP;
    else
        dwi_fov = DWIparams.FOV.RL;
        dwi_offc = DWIparams.Offc.RL;
    end
    dwi_th = DWIparams.slice_thickness;
    dwi_ovs_fac = nSLICE / DWIparams.nSLICE;
    dwi_fov_ovs = dwi_fov * dwi_ovs_fac;
    
    %* DWI slice center.
    %** Slice 1 to end is from +ve to -ve coord.
    %dwi_mm_v = linspace(-dwi_fov_ovs/2 + dwi_th/2,dwi_fov_ovs/2 - dwi_th/2,nSLICE);
    %dwi_mm_v = linspace(-dwi_fov/2 + dwi_th/2,dwi_fov/2 - dwi_th/2,nSLICE);
    dwi_mm_v = linspace(dwi_fov/2 - dwi_th/2,-dwi_fov/2 + dwi_th/2,nSLICE);
    dwi_mm_v = dwi_mm_v + dwi_offc; % DWI slice position in gradient coord   
    
    %* Check slice thickness.
    val = unique(abs(diff(dwi_mm_v)));
    if val~=dwi_th
        error('s_reslice_map_v2:main', ...
            'DWI slice thickness must be the same as slice center locations')
    end
    
    
    % Reslice SMAP volume at DWI position. This sl_smap_target_v is not
    % considered for the ovelap between chunks. It is considered below when
    % 3D multichunk acquisition is used.
    ref_slice_v = 1:snz;
    sl_smap_target_v = interp1(ref_mm_v, ref_slice_v, dwi_mm_v);
    if length(sl_smap_target_v)~=nSLICE
        error('s_reslice_map:main','sl_smap_target_v or sl_smap_ref_v must equal to nSLICE.')
    end
    
    %***
    % These two must be similar.
    %(sl_smap_target_v(end)-sl_smap_target_v(1))*ref_th
    %dwi_fov
    %***
        
    % Report.
    fprintf('      sl_smap_target_v\n')
    %disp(sl_smap_target_v)
    disp([sl_smap_target_v(1) sl_smap_target_v(end)])
    
    
    % Take care of 3D case where there is overlapping between chunks.
    % SMAP positions are repeated at overlapping slices. This may be
    % significant if there are many chunks.
    if strcmpi(DWIparams.scan_mode,'3D')
        nSlice = DWIparams.nSLICE;        
        nKz_per_chunk = RECONparams.ENCima.oversample_resolutions(3);
        nChunk = DWIparams.chunks;
        nSlice_keep = nSlice / nChunk;
        ovs_fac = RECONparams.ENCima.oversample_factors(3);
        
        nCutSlice_per_chunk = (nKz_per_chunk*nChunk - nSlice) / nChunk;
        % == (nSlice*ovs_fac - nSlice) / nChunk, cutoff slices per chunk as
        % much as each chunk is oversampled.
        
        nCutSlice_at_chunk_start = nCutSlice_per_chunk / 2;
        nCutSlice_at_chunk_start = floor(nCutSlice_at_chunk_start);
        nCutSlice_at_chunk_end = nKz_per_chunk - nCutSlice_at_chunk_start - nSlice_keep;
        nOverlapSlice_at_chunk_end = (nKz_per_chunk - nSlice_keep - ...
            nCutSlice_at_chunk_start) + (nCutSlice_at_chunk_start); % @(current chunk) + @(previous chunk)
        
        % Number of slices overlapped at the end of each chunk in
        % chunk #2 ~ nChunk. For example,
        %
        %        -----------------------------  chunk#: nChunk
        % Head   - 1 - 2 - 3 - 4 - 5 - 6 - 7 -   (slice index)              Foot
        %        -----------------------------
        %                        -----------------------------  chunk#: nChunk-1
        %                        - 1 - 2 - 3 - 4 - 5 - 6 - 7 -
        %                        -----------------------------
        %
        % This scheme is recognized from the study,
        % [scan20120703_r4651__Jeong_999999].
        
         
        % Consider redundant slices.
        % This is for [s_reslice_map_3d_v2.m] and 'Method 2-2' below.
        
        % SMAP slice position must be as below,
        %   - 1 - 2 - 3 - 4 - 5 - 6 - 7 -       DWI chunk n
        %                   - 1 - 2 - 3 - 4 - 5 - 6 - 7 -  DWI chunk n-1
        %
        %   - 1 - 2 - 3 - 4 - 5 - 6 - 7 -   SMAP for DWI chunk n
        %                   - 5 - 6 - 7 - 8 - 9 - 10 - 11 - SMAP for DWI chunk n-1
        %                                   - 9 - 10 - 11 - 12 - 13 - 14 - 14 - SMAP for DWI chunk n-2
        %
        % Then SMAP slices are
        %   1,2,3,4,5,6,7 | 5,6,7,8,9,10,11 | 9,10,11,12,13,14,15 | ...
        % then at each chunk position, there must be nSlice_keep continuous
        % SMAP slices and they must be also continuous slices between
        % chunks, such as | 2,3,4,5 | 6,7,8,9 | 10,11,12,13 | ... etc.
        %
        % The redundant SMAP slices per each chunk are as many as
        % nOverlapSlice_at_chunk_end. This is in 'Method 2-2' below.
        
        
        %----- Method 1 -----
%         index_sl_smap_target_overlap_v = [];
%         start_v = 1:(nKz_per_chunk-1):(nKz_per_chunk-1)*nChunk;
%         end_v = (nKz_per_chunk-1):(nKz_per_chunk-1):(nKz_per_chunk-1)*nChunk;
%         for ind = 1:length(start_v)
%             v = start_v(ind):end_v(ind);
%             v(end+1) = v(end);
%             index_sl_smap_target_overlap_v = [index_sl_smap_target_overlap_v, v];
%         end
%         sl_smap_target_v = sl_smap_target_v(index_sl_smap_target_overlap_v);
%         clear  index_sl_smap_target_overlap_v

        
        %----- Method 2 -----
%         index_sl_smap_target_overlap_v = [];
%         for ind_chunk = 1:nChunk
%             v = 1+nKz_per_chunk*(ind_chunk-1):nKz_per_chunk*ind_chunk;
%             v = v - nOverlapSlice_at_chunk_end*(ind_chunk-1);
%             index_sl_smap_target_overlap_v = [index_sl_smap_target_overlap_v,v];
%         end        
%         sl_smap_target_v = sl_smap_target_v(index_sl_smap_target_overlap_v);
%         clear  index_sl_smap_target_overlap_v
        

        %----- Method 2-2 -----
        index_sl_smap_target_overlap_v = [];
        v = zeros(1,nKz_per_chunk);
        for ind_chunk = 1:nChunk
            if ind_chunk==1
                v(1:nCutSlice_at_chunk_start) = (1-nCutSlice_at_chunk_start):0;
                v(nCutSlice_at_chunk_start+1:end) = 1:(nKz_per_chunk-nCutSlice_at_chunk_start);
            else
                v = 1+nKz_per_chunk*(ind_chunk-1):nKz_per_chunk*ind_chunk;
                v = v - (nOverlapSlice_at_chunk_end*(ind_chunk-1)+nCutSlice_at_chunk_start);
            end
            index_sl_smap_target_overlap_v = [index_sl_smap_target_overlap_v,v];
        end
        sl_smap_target_v = interp1(1:length(sl_smap_target_v), ...
            sl_smap_target_v,index_sl_smap_target_overlap_v,'linear','extrap');
        clear  index_sl_smap_target_overlap_v
                
        %*** BEWARE !!! ***
        % Actual slice coverage in DWI is dwi_fov_ovs / ovs_fac mm. Then
        %
        % (sl_smap_target_v(end)-sl_smap_target_v(1))*ref_th -
        % nCutSlice_at_chunk_start*dwi_th ~ dwi_fov_ovs / ovs_fac 
        %
        % all in mm unit.
        %
        % ******
                
        % Reserve the multi-chunk related parameters.
        DWIparams.multichunk.nCutSlice_per_chunk = nCutSlice_per_chunk;
        DWIparams.multichunk.nCutSlice_at_chunk_start = nCutSlice_at_chunk_start;
        DWIparams.multichunk.nCutSlice_at_chunk_end = nCutSlice_at_chunk_end;
        DWIparams.multichunk.nOverlapSlice_at_chunk_end = nOverlapSlice_at_chunk_end;
        DWIparams.multichunk.sl_smap_target_v = sl_smap_target_v;
        DWIparams.multichunk.string = ...
            ['nCutSlice_per_chunk: number of slices to be removed from   '; ...
             'each chunk, including front and end slices, before         '; ...
             'combining slices from all chunks.                          '; ...
             '                                                           '; ...
             'nCutSlice_at_chunk_start: number of slices to be removed   '; ...
             'from the start of each chunk before combining slices from  '; ...
             'all chunks.                                                '; ...
             '                                                           '; ...
             'nCutSlice_at_chunk_end: similarly at the end of each chunk '; ...
             '                                                           '; ...
             'nOverlapSlice_at_chunk_end: number of slices overlapped    '; ...
             'between the end of nth and start of (n-1)th chunks.        '; ...
             '                                                           '; ...
             'sl_smap_target_v: SMAP slice positions at DWI slice        '; ...
             'position considering slice overlap between chunks.         '];
    end
    
    
    % Save runtime parameters.
    % This is for updated DWIparams.
    cd(saveReconDir_s)
    save('params','DWIparams','-append')
    fprintf('      DWIparams.multichunk is updated\n')
    multichunk = DWIparams.multichunk;
    save DWIparams.multichunk.mat  multichunk
        
    
    % Report.
    fprintf('      sl_smap_target_v with 3D chunk overlapping\n')
    %disp(sl_smap_target_v)
    disp([sl_smap_target_v(1) sl_smap_target_v(end)])
    
    
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
    
    
    % Check if I_smap_img_4d can be saved as .mat file.
    filesize__ref = 2.0; % reference filesize in GB
    SMAPparams.filesize__ref = filesize__ref; 
    filesize__I_smap_img_4d = prod([sny,snx,length(sl_smap_target_v),snc])*8; % bytes
    SMAPparams.filesize__I_smap_img_4d = filesize__I_smap_img_4d/10^9; % in GB
    if SMAPparams.filesize__I_smap_img_4d > SMAPparams.filesize__ref
        flag_save_as_chunk = true;
    else
        flag_save_as_chunk = false;
    end
    
    
    % Save runtime SMAPparams.
    cd(saveReconDir_s)
    save('params','SMAPparams','-append')
    fprintf('      [SMAPparams.filesize__ref] and [SMAPparams.filesize__I_smap_img_4d] are updated\n')
    save  SMAPparams.filesize__ref.mat  filesize__ref
    save  SMAPparams.filesize__I_smap_img_4d.mat  filesize__I_smap_img_4d
    
    
    % Resample SMAP slices
    fprintf('    Resample SMAP slices\n')    
    [x,y,z] = meshgrid(1:snx,(1:sny)',1:snz);
    
    if flag_save_as_chunk==false
        
        I_smap_resample_4d = zeros(sny,snx,length(sl_smap_target_v),snc);
        for ind_coil = 1:snc
            smap_vol = I_smap_4d(:,:,:,ind_coil);
            [xi,yi,zi] = meshgrid(1:snx,(1:sny)',sl_smap_target_v);
            I_smap_resample_4d(:,:,:,ind_coil) = interp3(x,y,z,smap_vol,xi,yi,zi);
            if ind_coil==1
                I_mask_resample_3d = interp3(x,y,z,I_mask_S_3d,xi,yi,zi);
            end
        end
        
        % Remove holes in mask.
        for ind_slice = 1:size(I_mask_resample_3d,3)
            I_mask_resample_3d(:,:,ind_slice) = imfill(I_mask_resample_3d(:,:,ind_slice),'holes');
        end
        I_mask_resample_3d = (I_mask_resample_3d > 0);
        
        % Reserve resampled SMAP and MASK.
        %I_smap_nav_4d = I_smap_resample_4d;
        I_smap_img_4d = I_smap_resample_4d;
        %I_mask_S_nav_3d = I_mask_resample_3d;
        I_mask_S_img_3d = I_mask_resample_3d;
        clear  I_smap_resample_4d  I_mask_resample_3d  x  y  z  xi  yi  zi
        
        % Save SMAP.
        cd(saveDataDir_s)
        %eval(sprintf('save  I_smap_nav_4d    I_smap_nav_4d'))
        %eval(sprintf('save  I_mask_S_nav_3d  I_mask_S_nav_3d'))
        eval(sprintf('save  I_smap_img_4d    I_smap_img_4d'))
        eval(sprintf('save  I_mask_S_img_3d  I_mask_S_img_3d'))
        
    else % save as chunk
        
        I_smap_resample_4d = zeros(sny,snx,DWIparams.kz,snc);
        start_v = length(sl_smap_target_v)-DWIparams.kz+1:-DWIparams.kz:1;
        end_v = length(sl_smap_target_v):-DWIparams.kz:1;
        if length(start_v)~=DWIparams.chunks || length(end_v)~=DWIparams.chunks
            error('s_reslice_map_3d_v3:main','Check kz and number of chunks')
        end
        for ind_chunk = DWIparams.chunks:-1:1
            kz_v = start_v(ind_chunk):end_v(ind_chunk);
            sl_smap_target_chunk_v = sl_smap_target_v(kz_v);
            [xi,yi,zi] = meshgrid(1:snx,(1:sny)',sl_smap_target_chunk_v);
            for ind_coil = 1:snc
                smap_vol = I_smap_4d(:,:,:,ind_coil);
                I_smap_resample_4d(:,:,:,ind_coil) = interp3(x,y,z,smap_vol,xi,yi,zi);
                if ind_coil==1
                    I_mask_resample_3d = interp3(x,y,z,I_mask_S_3d,xi,yi,zi);
                    
                    % Remove holes in mask.
                    for ind_slice = 1:size(I_mask_resample_3d,3)
                        I_mask_resample_3d(:,:,ind_slice) = imfill(I_mask_resample_3d(:,:,ind_slice),'holes');
                    end
                    I_mask_resample_3d = (I_mask_resample_3d > 0);
                end
            end
            I_smap_img_4d = I_smap_resample_4d;
            I_mask_S_img_3d = I_mask_resample_3d;
            
            cd(saveDataDir_s)
            smapname_s = sprintf('I_smap_img_chunk%.2d_4d',ind_chunk);
            maskname_s = sprintf('I_mask_S_img_chunk%.2d_3d',ind_chunk);            
            eval(sprintf('save  %s  I_smap_img_4d',smapname_s))
            eval(sprintf('save  %s  I_mask_S_img_3d',maskname_s))
            fprintf('      %s and %s are saved\n',smapname_s,maskname_s)            
        end        
    end
    
    % Get new size.
    [sny,snx,snz,snc] = size(I_smap_img_4d);
    
        
    % Report.
    fprintf('\n')
    fprintf('      SMAP to be used are resampled and saved\n')
    if flag_save_as_chunk==false
        %fprintf('        I_smap_nav_4d,     [sny,snx,snz,snc]  = [%d,%d,%d,%d]\n',size(I_smap_nav_4d))
        %fprintf('        I_mask_S_nav_3d,   [sny,snx,snz]      = [%d,%d,%d]\n',size(I_mask_S_nav_3d))
        fprintf('        I_smap_img_4d,     [sny,snx,snz,snc]  = [%d,%d,%d,%d]\n',size(I_smap_img_4d))
        fprintf('        I_mask_S_img_3d,   [sny,snx,snz]      = [%d,%d,%d]\n',size(I_mask_S_img_3d))
    else
        %fprintf('        I_smap_nav_4d,     [sny,snx,snz,snc]  = [%d,%d,%d,%d]\n',size(I_smap_nav_4d))
        %fprintf('        I_mask_S_nav_3d,   [sny,snx,snz]      = [%d,%d,%d]\n',size(I_mask_S_nav_3d))
        fprintf('        I_smap_img_4d/chunk     [sny,snx,snz,snc]  = [%d,%d,%d,%d]\n',size(I_smap_img_4d))
        fprintf('        I_mask_S_img_3d/chunk   [sny,snx,snz]      = [%d,%d,%d]\n',size(I_mask_S_img_3d))
    end
    
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
    
    % Flip dim 3.
    I_fmap_temp_4d = I_fmap_4d*0;
    if strcmpi(DWIparams.scan_mode,'3D')
        for ind = 1:fnt
            I_fmap_temp_4d(:,:,:,ind) = flipdim(I_fmap_4d(:,:,:,ind),3);
        end
    end
    I_fmap_4d = I_fmap_temp_4d;
    clear  I_fmap_temp_4d
            
    
    
    % Get corresponding slice positions for FMAP
    fprintf('    Get corresponding FMAP slice positions\n')
    
    % Set corresponding slices.
    % Set corresponding DWI slices.
    if strcmpi(DWIparams.scan_mode,'3D')
        nSLICE = RECONparams.ENCima.oversample_resolutions(3) * DWIparams.chunks;
    else
        nSLICE = DWIparams.nSLICE;
    end
    nSLICE = DWIparams.nSLICE;
    
    
    % ---------- New routine for reslicing SMAP ----------
    
    % FMAP slice FOV and Offc.
    %* Assume FMAP data is in the same FOV and Offc as DWI data.
    if strcmpi(REFparams.slice_orientation,'transverse')
        fmap_fov = FMAPparams.FOV.FH;
        fmap_offc = DWIparams.Offc.FH;
    elseif strcmpi(REFparams.slice_orientation,'coronal')
        fmap_fov = FMAPparams.FOV.AP;
        fmap_offc = DWIparams.Offc.AP;
    else
        fmap_fov = FMAPparams.FOV.RL;
        fmap_offc = DWIparams.Offc.RL;
    end
    fmap_th = FMAPparams.slice_thickness;
    fmap_ovs_fac = 1; % 1 since FMAP is from PAR/REC file
    fmap_fov_ovs = fmap_fov * fmap_ovs_fac;
    
    %* REF slice center.
    %fmap_mm_v = linspace(-fmap_fov_ovs/2 + fmap_th/2,fmap_fov_ovs/2 - fmap_th/2,fnz);
    fmap_mm_v = linspace(fmap_fov_ovs/2 - fmap_th/2,-fmap_fov_ovs/2 + fmap_th/2,fnz);
    fmap_mm_v = fmap_mm_v + fmap_offc; % FMAP slice position in grad coord
    
    %* Check slice thickness.
    val = unique(abs(diff(fmap_mm_v)));
    if val~=fmap_th
        error('s_reslice_map_v2:main', ...
            'FMAP slice thickness must be the same as slice center locations')
    end
    
    
    % Reslice FMAP volume at DWI position. This sl_fmap_target_v is not
    % considered for the ovelap between chunks. It is considered below when
    % 3D multichunk acquisition is used.
    fmap_slice_v = 1:fnz;
    sl_fmap_target_v = interp1(fmap_mm_v, fmap_slice_v, dwi_mm_v);
    if length(sl_fmap_target_v)~=nSLICE
        error('s_reslice_map:main','sl_fmap_target_v or sl_fmap_ref_v must equal to nSLICE.')
    end
    
    %***
    % These two must be similar.
    %(sl_fmap_target_v(end)-sl_fmap_target_v(1))*fmap_th
    %dwi_fov
    %***
    
    
    % Report.
    fprintf('      sl_fmap_target_v\n')
    disp([sl_fmap_target_v(1) sl_fmap_target_v(end)])
    
    
    % Take care of 3D case where there is overlapping between chunks.
    % FMAP positions are repeated at overlapping slices. This may be
    % significant if there are many chunks.
    if strcmpi(DWIparams.scan_mode,'3D')
        nSlice = DWIparams.nSLICE;
        nKz_per_chunk = RECONparams.ENCima.oversample_resolutions(3);
        nChunk = DWIparams.chunks;
        nSlice_keep = nSlice / nChunk;
        ovs_fac = RECONparams.ENCima.oversample_factors(3);
        
        nCutSlice_per_chunk = (nKz_per_chunk*nChunk - nSlice) / nChunk;
        % == (nSlice*ovs_fac - nSlice) / nChunk, cutoff slices per chunk as
        % much as each chunk is oversampled.
        
        nCutSlice_at_chunk_start = nCutSlice_per_chunk / 2;
        nCutSlice_at_chunk_start = floor(nCutSlice_at_chunk_start);
        nCutSlice_at_chunk_end = nKz_per_chunk - nCutSlice_at_chunk_start - nSlice_keep;
        nOverlapSlice_at_chunk_end = (nKz_per_chunk - nSlice_keep - ...
            nCutSlice_at_chunk_start) + (nCutSlice_at_chunk_start); % @(current chunk) + @(previous chunk)
        
        % Number of slices overlapped at the end of each chunk in
        % chunk #2 ~ nChunk. For example,
        %
        %        -----------------------------  chunk#: nChunk
        % Head   - 1 - 2 - 3 - 4 - 5 - 6 - 7 -   (slice index)              Foot
        %        -----------------------------
        %                        -----------------------------  chunk#: nChunk-1
        %                        - 1 - 2 - 3 - 4 - 5 - 6 - 7 -
        %                        -----------------------------
        %
        % This scheme is recognized from the study,
        % [scan20120703_r4651__Jeong_999999].
        
        
        % Consider redundant slices.
        % This is for [s_reslice_map_3d_v2.m] and 'Method 2-2' below.
        
        % SMAP slice position must be as below,
        %   - 1 - 2 - 3 - 4 - 5 - 6 - 7 -       DWI chunk n
        %                   - 1 - 2 - 3 - 4 - 5 - 6 - 7 -  DWI chunk n-1
        %
        %   - 1 - 2 - 3 - 4 - 5 - 6 - 7 -   FMAP for DWI chunk n
        %                   - 5 - 6 - 7 - 8 - 9 - 10 - 11 - FMAP for DWI chunk n-1
        %                                   - 9 - 10 - 11 - 12 - 13 - 14 - 14 - FMAP for DWI chunk n-2
        %
        % Then FMAP slices are
        %   1,2,3,4,5,6,7 | 5,6,7,8,9,10,11 | 9,10,11,12,13,14,15 | ...
        % then at each chunk position, there must be nSlice_keep continuous
        % SMAP slices and they must be also continuous slices between
        % chunks, such as | 2,3,4,5 | 6,7,8,9 | 10,11,12,13 | ... etc.
        %
        % The redundant SMAP slices per each chunk are as many as
        % nOverlapSlice_at_chunk_end. This is in 'Method 2' below.
        
        
        %----- Method 1 -----
        %         index_sl_smap_target_overlap_v = [];
        %         start_v = 1:(nKz_per_chunk-1):(nKz_per_chunk-1)*nChunk;
        %         end_v = (nKz_per_chunk-1):(nKz_per_chunk-1):(nKz_per_chunk-1)*nChunk;
        %         for ind = 1:length(start_v)
        %             v = start_v(ind):end_v(ind);
        %             v(end+1) = v(end);
        %             index_sl_smap_target_overlap_v = [index_sl_smap_target_overlap_v, v];
        %         end
        %         sl_smap_target_v = sl_smap_target_v(index_sl_smap_target_overlap_v);
        %         clear  index_sl_smap_target_overlap_v
        
        
        %----- Method 2 -----
        %         index_sl_smap_target_overlap_v = [];
        %         for ind_chunk = 1:nChunk
        %             v = 1+nKz_per_chunk*(ind_chunk-1):nKz_per_chunk*ind_chunk;
        %             v = v - nOverlapSlice_at_chunk_end*(ind_chunk-1);
        %             index_sl_smap_target_overlap_v = [index_sl_smap_target_overlap_v,v];
        %         end
        %         sl_smap_target_v = sl_smap_target_v(index_sl_smap_target_overlap_v);
        %         clear  index_sl_smap_target_overlap_v
        
        
        %----- Method 2-2 -----
        index_sl_fmap_target_overlap_v = [];
        v = zeros(1,nKz_per_chunk);
        for ind_chunk = 1:nChunk
            if ind_chunk==1
                v(1:nCutSlice_at_chunk_start) = (1-nCutSlice_at_chunk_start):0;
                v(nCutSlice_at_chunk_start+1:end) = 1:(nKz_per_chunk-nCutSlice_at_chunk_start);
            else
                v = 1+nKz_per_chunk*(ind_chunk-1):nKz_per_chunk*ind_chunk;
                v = v - (nOverlapSlice_at_chunk_end*(ind_chunk-1)+nCutSlice_at_chunk_start);
            end
            index_sl_fmap_target_overlap_v = [index_sl_fmap_target_overlap_v,v];
        end
        sl_fmap_target_v = interp1(1:length(sl_fmap_target_v), ...
            sl_fmap_target_v,index_sl_fmap_target_overlap_v,'linear','extrap');
        clear  index_sl_fmap_target_overlap_v
        
        
        %*** BEWARE !!! ***
        % Actual slice coverage in DWI is dwi_fov_ovs / ovs_fac mm. Then
        %
        % (sl_fmap_target_v(end)-sl_fmap_target_v(1))*fmap_th - nCutSlice_at_chunk_start*dwi_th ~ 
        % dwi_fov_ovs / ovs_fac
        %
        % all in mm unit.
        %
        % ******
        
        
        % No need for reserving the multi-chunk related parameters.
%         DWIparams.multichunk.nCutSlice_per_chunk = nCutSlice_per_chunk;
%         DWIparams.multichunk.nCutSlice_at_chunk_start = nCutSlice_at_chunk_start;
%         DWIparams.multichunk.nCutSlice_at_chunk_end = nCutSlice_at_chunk_end;
%         DWIparams.multichunk.nOverlapSlice_at_chunk_end = nOverlapSlice_at_chunk_end;
%         DWIparams.multichunk.sl_smap_target_v = sl_smap_target_v;
%         DWIparams.multichunk.string = ...
%             ['nCutSlice_per_chunk: number of slices to be removed from   '; ...
%             'each chunk, including front and end slices, before         '; ...
%             'combining slices from all chunks.                          '; ...
%             '                                                           '; ...
%             'nCutSlice_at_chunk_start: number of slices to be removed   '; ...
%             'from the start of each chunk before combining slices from  '; ...
%             'all chunks.                                                '; ...
%             '                                                           '; ...
%             'nCutSlice_at_chunk_end: similarly at the end of each chunk '; ...
%             '                                                           '; ...
%             'nOverlapSlice_at_chunk_end: number of slices overlapped    '; ...
%             'between the end of nth and start of (n-1)th chunks.        '; ...
%             '                                                           '; ...
%             'sl_smap_target_v: SMAP slice positions at DWI slice        '; ...
%             'position considering slice overlap between chunks.         '];
    end

    
    % Report.
    fprintf('      sl_fmap_target_v with 3D chunk overlapping\n')
    %disp(sl_fmap_target_v)
    disp([sl_fmap_target_v(1) sl_fmap_target_v(end)])
    
    
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
clear  I_*

% Pack memory.
pack

fprintf('\n\n')










