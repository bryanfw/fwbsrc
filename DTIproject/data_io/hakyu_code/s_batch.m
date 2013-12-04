
%[s_batch.m]
%
%
% Ha-Kyu



%% Set directory
clear all
close all
pack
clc
warning off

disp('-------------------- Set directories --------------------')

% Study parameters.
studyname_s = 'scan20120723_r4767__Anderson_307740';
examcardFilename_s = 'ExamCard_to_XML_201207231657014535000.xml';
anatomy_c = {'optic_nerve','spinal_cord','brain','leg','phantom','auditory_canal'};
anatomy_s = anatomy_c{3};

% Define directory.
if strcmpi(computer,'PCWIN') || strcmpi(computer,'PCWIN64')
    projectDir_s = 'Y:\Scan';
    %projectDir_s = 'X:\3T\ON\controls';
    scriptDir_s = sprintf('Z:\\Research\\PPE\\DTI_NAV\\Matlab\\%s\\%s\\',...
        anatomy_s,studyname_s);
else
    % Add .m file paths.
    cd('/home/jeongh1')
    hkj_addpath

    projectDir_s = '/home/DTI_NAV/Scan';
    %projectDir_s = '/home/ms_sas/3T/ON/controls';
    scriptDir_s = sprintf('/home/jeongh1/Research/PPE/DTI_NAV/Matlab/%s/%s/', ...
        anatomy_s, studyname_s);
end

loadDir_s = sprintf('%s%s%s%slistdata%s', ...
    projectDir_s,filesep,studyname_s,filesep,filesep);
sharedDir_s = sprintf('%s%s%s%slistdata%sshared%s', ...
    projectDir_s,filesep,studyname_s,filesep,filesep,filesep);
loadParDir_s = sprintf('%s%s%s%sparrec%s', ...
    projectDir_s,filesep,studyname_s,filesep,filesep);
examcardDir_s = sprintf('%s%s%s%sexamcard%s', ...
    projectDir_s,filesep,studyname_s,filesep,filesep);

% Move to script directory.
cd(scriptDir_s)
fprintf('\n\n')



%% Run
disp('-------------------- Run procedures --------------------')
for ind_raw_file = [6,7] % [6,7] raw DWI file index as shown in the raw file name
    
    
    
    %-------------------- Batch process --------------------
    % Control overal batch process
    BATCHparams.read_data = true;
    BATCHparams.recon_data = true;
    
    % Get SSH PAR/REC filename.
%     SSHparrecfilename = {'Anderson_306670_WIP_sshDTI_SENSE_13_1.PAR',...
%         'Anderson_306670_WIP_sshDTI_SENSE_14_1.PAR',...
%     'Anderson_306670_WIP_sshDTI_SENSE_15_1.PAR'}; % or {}
    SSHparrecfilename = {};
    
    % Get FMAP PAR/REC filename.
    %* 'Smith_4434_WIP_B0_Map_12_1.PAR' has problem, not a valid FMAP.
    %* Check if 'B0 field map' of 'calibrate' can cause this problem.
    FMAPparrecfilename = ['70322_Anderson_307440_WIP_FMAP_CLEAR_5_1.PAR']; % or []
    
        
    % Get REF and DWI raw filename and scan index.
    % When there are scans stopped during study, the scan index in the
    % ExamCard may be different than actual scan number. So enter both of
    % the raw data file name and examcard index.
    REFrawfilename = 'raw_004';
    REFexamcardindex = 4; % corresponding index for the scan in the ExamCard, e.g., 7th scan in ExamCard
    DWIrawfilename = sprintf('raw_%.3d',ind_raw_file); 
    if ind_raw_file==6
        DWIexamcardindex = ind_raw_file;
    else
        DWIexamcardindex = ind_raw_file+1;
    end
    
    
        
    %-------------------- Scanner parameter --------------------
    % Main field strength.
    main_field_strength = 70000; % Gauss
    
    % Coil type: 'Coil selection' parameter in GEOMETRY tab.
    % This is for using 32 channel coils in future. With 32 channel coil,
    % the parameter orders in the ExamCard will change.
    
    % Typically,
    % 7T:   'SENSE-Head-7T' - 16 channel
    %       'RX-Intf-1' - 32 channel
    %       'RX-Intf-1_Quad-TR-1' - 16 channel spine coil
    % 3T:   'SENSE-Head-8' - 8 channel
    %       'SENSE-NV-16' - 16 channel for spinal cord
    %       'SENSE-Breast-4' - 4 channel for breast
    %       'SENSE-Head-32' - 32 channel
    EX_GEO_cur_stack_coil_id = 'RX-Intf-1';
    coilType = EX_GEO_cur_stack_coil_id;
     
    
    
    
    %-------------------- Patient geometry --------------------
    % Patient orientation: supine, prone, rightdecubitus, leftdecubitus    
    patient_orientation = 'supine';
    
    % Patient position: headfirst, feetfirst
    patient_position = 'headfirst';
    
    
    
        
    %-------------------- Data processing --------------------
    
    %===== DWIparams =====
    
    % Define user_defined_dw_ori. This goes into [f_get_dw_orientation3.m].
    % 
    % User defined order entered in 'scanner console' is [L,P,H] order. But
    % that needs to be changed from LPH system to PHL system as shown in
    % PAR file before it is feeded into [f_get_dw_orientation3.m] within
    % [s_get_DWI_RECONparam.m].
    
    switch ind_raw_file       
        otherwise
            user_defined_dw_ori = []; % or [] in [RL,AP,FH], dim = [ndir,3]
            if ~isempty(user_defined_dw_ori)
                user_defined_dw_ori = user_defined_dw_ori(:,[2,3,1]); % from [LPH] to [PHL]
            end
    end
    
    
    
    %===== SMAPparams =====
    
    % This parameter is for coil sensitivity generation.
    smap_fit_method_c = {'l2p','tps'}; %[local 2D polynomial fit, thin-plate spline fit]
    smap_fit_method = smap_fit_method_c{2};
    
    % Generate bodycoil equivalent using reference coil images.
    %* When bodycoil data are not available or reference data are to be
    %* used for bodycoil equivalent, modify this parameter.
    %* 'body': use bodycoil (or volumecoil) data as bodycoil equivalent
    %* 'ref':  use sum-of-square of reference coil data for bodycoil
    %*         equivalent.
    smap_gen_bodycoil = 'body'; % ['body','ref'] 
   
    
    
    %===== DATAFLAGparams =====
    
    % This parameter determines reading and saving raw files and generating
    % SMAP and FMAP.
    
    % Use total variation in selecting the best EPI ghost correction
    % results.
    %* This is implemeted in [f_epi_phase_corr_phc_v2.m] in
    %* [s_read_DWI6.m].
    %* It work for phantom where SNR is pretty high. But may not work for
    %* in vivo acquisitions. 
    epi_corr_sel_total_var = 0; %[0,1] = [don't, do] --- always 0
    
    % Set EPI echo phase correction method.
    %* If 'epi_corr_sel_total_var' is 0, then this will be used.
    %* If epi_corr_fit_method == 0 then no correction is done.
    switch coilType
        case 'SENSE-Breast-4'
            epi_corr_fit_method = 4;
        case 'RX-Intf-1_Quad-TR-1'
            epi_corr_fit_method = 4; % or use 0, don't use 3
            
            %--- Note ---
            % For [scan20120420_r4140__Anderson_306871], see EPI phase
            % correction results for slice 6-coil 9-avg 1: Method 4 is much
            % better than all the others, similar as no correction. slice
            % 6-coil 11-avg 1: Method 4 is slightly better than 2, much
            % better than 1 and 3, slightly better than no correction.
            
        case 'RX-Intf-1'
            %epi_corr_fit_method = 1; % linear fit
            epi_corr_fit_method = 1;
        case 'SENSE-Head-7T'
            epi_corr_fit_method = 1;
        case 'SENSE-Head-8'
            epi_corr_fit_method = 1;
        case 'SENSE-NV-16'
            epi_corr_fit_method = 1;
        case 'SENSE-Head-32'
            epi_corr_fit_method = 1;
        otherwise
            error('Unknown coilType')
    end
    
    % Match object position between SMAP and IMG data.
    %* Use centroids difference between reconstructed then masked IMG image
    %* and masked SMAP image.
    match_obj_pos_smap_img = 1; %[0,1]=[don't,do] --- preferably 1
    
    % Set raw data related process.
    %calc_ssh_diff = false;
    switch ind_raw_file
        case 6
            read_ref = ~BATCHparams.read_data; % read and save coil reference data (REF)
            read_dwi = ~BATCHparams.read_data; % read and save dwi data (DWI)
            gen_fmap = ~BATCHparams.read_data; % generate B0 fieldmap (FMAP)
            gen_smap = ~BATCHparams.read_data; % generate coil sensitivity map (SMAP)
        case 7
            read_ref = ~BATCHparams.read_data;
            read_dwi = ~BATCHparams.read_data;
            gen_fmap = ~BATCHparams.read_data;
            gen_smap = ~BATCHparams.read_data;       
        otherwise
            error('Unknown ind_raw_file for DATAFLAGparams')
    end
    
    % Additional term at the end of, e.g., 'raw_011_R3E9S8' as
    % 'raw_011_R3E9S8_recon_again' etc.
    add_term_at_data_dir = [];% or []
    
    % NAV coregistration to IMG in [s_fmap_corr_NAV.m].
    coreg_nav_to_img = false; % should always be false
    
    
    
    
    %===== RECONFLAGparams =====
    
    % nav_fmap_corr and nav_deform are mutually exclusive except when both
    % are zero.
    phase_corr = true;      % DWI reconstruction with navigator phase correction
    nav_fmap_corr = true;
    if isempty(FMAPparrecfilename)
        nav_fmap_corr = false;
    end
    nav_deform = false;     % use deformation of nav to img, zero currently
    psi_method = 1;         % noise correlation method, default:1
    calc_diff = true;       % calculate diffusion properties, true or false
    
    % Use precalculation of encoding matrix.
    precalc_encoding_matrix_flag = false; % --- mostly false
    
    % Use sub-encoding matrix to expedite inverting the encoding matrix.
    % [0,1,2,3]=[Ny,R,Ns,R*Ns] index of column size of sub-encoding
    % matrix, Ny==R*Ns*EPI_FACTOR. This is sub_encoding_matrix_siz_idx.
    
    %*** IMPORTANT ***
    % When precalculation of encoding matrix is used
    % sub_encoding_matrix_siz_idx is better to be 0 (==Ny, whole column
    % inversion).
    %
    % It's been faster not using precalculation and sub-encoding matrix.
    %*** IMPORTANT ***
    
    sub_encoding_matrix_flag = false; % use sub-enc or not, --- mostly false
    sub_encoding_matrix_siz_idx = 0; % [0,1,2,3]==[Ny,R,Ns,R*Ns]
    if precalc_encoding_matrix_flag==true
        sub_encoding_matrix_flag = false;
        sub_encoding_matrix_siz_idx = 0;
    end
    
    % Use Parallel Computing Toolbox (PCT).
    %* Check if there is PCT.
    flag_pct = 0; % better not using PCT due to memory burden
    ind_pct = [];
    A = ver;
    for ind=1:length(A)
        s = A(ind).Name;
        if strcmpi(s,'Parallel Computing Toolbox')
            flag_pct = 1; % there is PCT
            ind_pct = ind;
            break
        end
    end
    
    %* Determine if PCT will be used.
    
    %*** IMPORTANT ***
    % Run PCT only when 'precalc_encoding_matrix_flag==true' is used.
    % Otherwise it will take longer time.
    %*** IMPORTANT ***
    
    run_pct = false; % will run PCT
    if flag_pct==1 && precalc_encoding_matrix_flag==true && run_pct==true
        run_pct = true;
    else
        run_pct = false;
    end
    
    % Resize to recon resolution or not.
    flag_resize_DWI = 1;
    
    % Average index: This is to recon image (in the final recon
    % [s_recon_DWI2.mat]) with single acquisition.
    flag_average_idx = []; % 1 ~ DWIparams.nNSA, if [] then use all avgs.
    
    
    
    
    %===== POSTparams =====
    
    % Co-registration.
    optim_s = 'powell'; % 'powell' (default) or 'simplex'
    method_s = 'rigid'; % 'rigid' or 'affine' or 'none'
    cost_s = 'nmi'; % 'nmi' (default),'mi','ncc','ecc' or 'cor'
    
    
    % Tensor calculation.
    medfilt_s = 'off'; % median filtering: 'on' or 'off'
    medfilt_kern_siz = 2; % size of [nxn] neighbor
    
    
    % ZOOM imaging.
    zoom_s = false; % [true,false],ZOOM imaging (IVI or OVS) is on
            
    

    
    %-------------------- Main procedures --------------------
    % Run process.
    t0_batch = clock;
    
    cd(scriptDir_s);
    s_recon_main
    
    t1_batch = etime(clock,t0_batch);
    fprintf('BATCH PROCESS FINISHED in %f min\n',t1_batch/60)
    
    if run_pct==1
        matlabpool close
        fprintf('\n')
        fprintf('PARALLEL COMPUTING TOOLBOX closes with ''local %d''\n',nworkers)
        fprintf('\n')
    end
    
end
clear all
warning on
pack




%% END










