
%[s_recon_main]
%
%
% Ha-Kyu



%% ----------------------------------------------------------------------------
%% Get ExamCard parameters
%   Get ExamCard parameters for use in,
%       REFparams
%       DWIparams
%% ----------------------------------------------------------------------------

%% Get ExamCard parameters

disp('-------------------- Get ExamCard parameters --------------------')

% Get ExamCard in XML format.
xmlfile = fullfile(examcardDir_s,examcardFilename_s);

% Read ExamCard parameter.
if main_field_strength==70000
    [refparams,dwiparams] = f_read_examcard_params_7T( ...
        xmlfile, REFrawfilename, DWIrawfilename, ...
        REFexamcardindex, DWIexamcardindex, coilType);
elseif main_field_strength==30000
    [refparams,dwiparams] = f_read_examcard_params_3T( ...
        xmlfile, REFrawfilename, DWIrawfilename, ...
        REFexamcardindex, DWIexamcardindex, coilType);
else
    error('Unknown main_field_strength')
end












%% ----------------------------------------------------------------------------
%% Get various parameters
%   GENparams - General
%   REFparams - Coil reference scan
%   DWIparams - DWI scan
%   RECONparams - DWI reconstruction
%   SMAPparams - Coil sensitivity map
%   FMAPparams - B0 fieldmap
%   DATAFLAGparams - Read REF, DWI, FMAP data or generate SMAP
%   RECONFLAGparams - Recon prarameters, e.g., phase correction, filedmap correction etc.
%   POSTparams - Post-processing, e.g., co-registration, tensor calculation etc.
%% ----------------------------------------------------------------------------


%% Get parameters
% Get paramters for scanning, reconstruction and post-processing.

%---------------------------------------------------------
% GENERAL scan parameters
disp('-------------------- GENparams --------------------')
s_get_GENparams
disp(GENparams)


%---------------------------------------------------------
% COIL REFERENCE scan parameters
disp('-------------------- REFparams --------------------')
s_get_REFparams
disp(REFparams)


%---------------------------------------------------------
% DWI and RECON parameters
disp('-------------------- DWIparams --------------------')
s_get_DWI_RECONparams
disp(DWIparams)

disp('-------------------- RECONparams --------------------')
disp(RECONparams)

disp('-------------------- Display INFO PAGE --------------------')
disp(RECONparams.ACQREC)


%---------------------------------------------------------
% Post-processing parameters
disp('-------------------- POSTparams --------------------')
s_get_POSTparams
disp(POSTparams)



disp('-------------------- Check SMAP and ALIASED matrix size --------------------')
s_check_smapSizeInfo

fprintf('\n\n')



%% Set coil sensitivity map generation related parameters
disp('-------------------- SMAPparams --------------------')
s_get_SMAPparams
disp(SMAPparams)



%% Set fieldmap postprocessing related parameters
disp('-------------------- FMAPparams --------------------')
s_get_FMAPparams
disp(FMAPparams)



%% Set up data processing related parameters
disp('-------------------- DATAFLAGparams --------------------')
s_get_DATAFLAGparams
disp(DATAFLAGparams)



%% Initialize recon flag parameters for final IMG recon
disp('-------------------- RECONFLAGparams --------------------')
s_get_RECONFLAGparams
disp(RECONFLAGparams)



%% Generate directories for recon data
disp('-------------------- Generate directories --------------------')
s_gen_dirs



%% Display POSTparams
disp('-------------------- POSTparams --------------------')
disp(POSTparams.coreg)
disp(POSTparams.tensor)



%% Other
disp('-------------------- Other params --------------------')
disp(flag_average_idx)



%% Clear unused data
clear  main_field_strength
clear  refparams
clear  dwiparams  DWIrawfilename  SSHparrecfilename  data  info
clear  epi_corr_fit_method  read_ref  read_dwi  gen_smap  gen_fmap
clear  phase_corr  nav_fmap_corr  nav_deform  psi_method  calc_diff
clear  patient_orientation  patient_position
fprintf('\n\n')



%% Save all parameters
cd(saveReconDir_s)
save('params','BATCHparams','DATAFLAGparams','DWIparams','FMAPparams',...
    'GENparams','POSTparams','RECONFLAGparams','RECONparams','REFparams',...
    'SMAPparams','examcardDir_s','loadDir_s','loadParDir_s','projectDir_s',...
    'saveDataDir_s','saveReconDir_s','scriptDir_s','sharedDir_s')
fprintf('*** params.mat is saved in [%s]\n',saveReconDir_s)

% Load runtime variables for updating param structs.
if exist('DATAFLAGparams.matchObjPosDelta.mat','file')
    load  DATAFLAGparams.matchObjPosDelta.mat
    DATAFLAGparams.matchObjPosDelta = matchObjPosDelta;
    fprintf('       DATAFLAGparams.matchObjPosDelta is loaded\n')
end
if exist('DWIparams.multichunk.mat','file')
    load  DWIparams.multichunk.mat
    DWIparams.multichunk = multichunk;
    fprintf('       DWIparams.multichunk is loaded\n')
end
if exist('SMAPparams.filesize__ref.mat','file')
    load  SMAPparams.filesize__ref.mat
    load  SMAPparams.filesize__I_smap_img_4d.mat
    SMAPparams.filesize__ref = filesize__ref;
    SMAPparams.filesize__I_smap_img_4d = filesize__I_smap_img_4d;
    fprintf('       SMAPparams.filesize_ref and .filesize__I_smap_img_4d are loaded\n')
end

fprintf('\n')





%% ----------------------------------------------------------------------------
%% Read LIST and PAR/REC data and process SMAP
%   REF
%   DWI
%   FMAP
%   SMAP
%% ----------------------------------------------------------------------------


%% Read, save and show REFERENCE and BODYCOIL data
disp('-------------------- Read REF --------------------')
s_read_REF



%% Generate coil sensitivity map in original REF size
disp('-------------------- Generate coil sensitivity map --------------------')
s_gen_SMAP



%% Read and resize fieldmap
disp('-------------------- Read FMAP --------------------')
s_read_FMAP



%% Read and show SENSE data
disp('-------------------- Read DWI --------------------')
s_read_DWI6 % takes care of scan mode = 3D



%% Clear unused data

clear  K_body  K_body_3d  K_nav*  K_img*  K_alias*  K_ref  K_ref_4d  
clear  I_smap  I_smap_3d  K_alias_3d  K_nav_3d
pack
fprintf('\n\n')






%% ----------------------------------------------------------------------------
%% Prepare SMAP and FMAP
%   This includes,
%       1. Adjust size of SMAP for reconstruction.
%       2. Generate in-registration state between S-F-NAV, so that NAV can be 
%          easily deformed into IMG using FMAP.
%       3. If necessary, SMAP can be deformed separately into NAV and IMG for 
%          more accurate reconstruction.
%% ----------------------------------------------------------------------------


%% Resize SMAP
disp('-------------------- Resize SMAP --------------------')
s_resize_SMAP



%% Reslice SMAP and FMAP slices for NAV and IMG recon
disp('-------------------- Reslice SMAP and FMAP --------------------')
if ~strcmpi(DWIparams.scan_mode,'3D')
    s_reslice_map
else
    s_reslice_map_3d_v3    
end






%% ----------------------------------------------------------------------------
%% Recon IMG and NAV
%   This includes,
%       1. Reconstruct IMG with no phase correction for comparison
%       2. Reconstruct NAV with no phase correction (of course) with and
%          without k-space window
%% ----------------------------------------------------------------------------


%% Recon and save IMG data for b-values and slices
disp('-------------------- Recon and save IMG data for b-values and slices --------------------')

% Check if IMG is already reconstructed.
cd(saveReconDir_s)
if ~( exist(sprintf('recon_img_ori00__R%dS%d.mat',R,Ns),'file')==2 ) && ...
        (BATCHparams.recon_data==1)
    
    % Recon IMG data for all slices.
    s_recon_IMG3    
end
fprintf('\n\n')



%% Recon and save NAV data for b-values, shots and slices
disp('-------------------- Recon and save NAV data for b-values, shots and slices --------------------')

% Check if NAV is already reconstructed.
cd(saveReconDir_s)
if ~( exist(sprintf('recon_nav_ori00__R%dS%d.mat',R,Ns),'file')==2 ) && ...
        ~( exist(sprintf('recon_nav_win_ori00__R%dS%d.mat',R,Ns),'file')==2 ) && ...
        (BATCHparams.recon_data==1)
    
    % Reconstruction of NAV.
    s_recon_NAV
end
fprintf('\n\n')



%% ----------------------------------------------------------------------------
%% Further processes on NAV
%   This may include,
%       1. Fieldmap correction of NAV
%       2. Direct deformation of NAV to IMG
%   NAV fieldmap correction and deformation will be combined together.
%% ----------------------------------------------------------------------------


%% Fieldmap correction of NAV for all nDW_GRAD and nSLICE
disp('-------------------- Fieldmap correction of NAV --------------------')

% Check if NAV is already FMAP corrected.
cd(saveReconDir_s)
if ~( exist(sprintf('recon_nav_win_fmap_ori00__R%dS%d.mat',R,Ns),'file')==2 ) && ...
        (~isempty(FMAPparams.filename))
    
    % Fieldmap correction of NAV.    
    s_fmap_corr_NAV_v2    
end
fprintf('\n\n')



%% ----------------------------------------------------------------------------
%% Reconstruct DWI IMG
%   This is for all b-values and slices with phase correction.
%% ----------------------------------------------------------------------------


%% SENSE recon DWI
disp('-------------------- Reconstruct DWI IMG with phase correction --------------------')

s_recon_DWI2



%% ----------------------------------------------------------------------------
%% Final resize of DWI IMG to recon voxel size
%   SENSE recon DWI images are resized to have the final recon voxel size
%% ----------------------------------------------------------------------------


%% Final resize of DWI recon images to recon voxel size
disp('-------------------- Final resize of DWI recon images to recon voxel size --------------------')
s_resize_DWI



%% END ALL PROCESSES
      



















