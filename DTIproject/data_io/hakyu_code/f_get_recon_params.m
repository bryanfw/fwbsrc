
function RECONparams = f_get_recon_params(REFparams,DWIparams,STACKima)
%[f_get_recon_params] calculates reconstruction related parameters such as 
% reconstruction matrix size and voxel size.
%
% USAGE
%   RECONparams = get_recon_params(REFparams,DWIparams,STACKima)
%
%
% Modified
% 2009.11.12.
% 2010.02.18. Add RECONparams.acq_vox_m  and RECONparams.nav__acq_vox_m.
%             Add varargin as verbose.
% 2010.07.09. [get_recon_params_v2.m].
%             Take STACK`ima.aspect_ratio as in input to [get_ENCima_v3.m].
% 2010.07.09. [get_recon_params_v3.m].
%             Recalculate the cell 'Calculate recon matrix and voxel size' for
%             rectangular FOV case also.
% 2010.07.12. [get_recon_params_v3.m].
%             Take REFparams as input.
% 2010.07.22. Modify to 'a_N_recon_p = double(uint8(E * ns * R));'.
%             This problem is found for processing foldover=AP,fatshift=P data.
% 2010.08.06.
%   Function name changed to from [get_recon_params_v3.m] to
%   [f_get_recon_params.m] with prefix 'f' to show that this is function.
% 2010.08.13.
%   Use DWIparams generated using [f_read_examcard_params.m].
% 2010.08.17.
%   Change R from uint8 to double.
% 2010.09.16.
%   Adjust final matrix size based on slice_orientation of REFparams and
%   DWIparams.
% 2010.09.24.
%   Remove varargin input argument.
% 2010.10.22.
%   Add DWIparams.shot_mode check.
%
% Ha-Kyu



%% Parse input
ACQ_RESOL = DWIparams.ACQ_RESOL;
RECON_RESOL = DWIparams.RECON_RESOL;
SENSE_FACTOR = DWIparams.SENSE_FACTOR;
isEPI = DWIparams.isEPI;
EPI_FACTOR1 = DWIparams.EPI_FACTOR(1);
if strcmpi(DWIparams.shot_mode,'multishot')
    EPI_FACTOR2 = DWIparams.EPI_FACTOR(2);
end
slice_thickness = DWIparams.slice_thickness;
if strcmpi(DWIparams.slice_orientation,'transverse')    
    if strcmpi(DWIparams.fold_over_dir,'AP')
        FOV = [DWIparams.FOV.RL, DWIparams.FOV.AP]; % [M,P]
    elseif strcmpi(DWIparams.fold_over_dir,'RL')
        FOV = [DWIparams.FOV.AP, DWIparams.FOV.RL]; % [M,P]
    else
        error('f_get_recon_params:main','Unknown DWI fold-over direction')
    end
elseif strcmpi(DWIparams.slice_orientation,'sagittal')    
    if strcmpi(DWIparams.fold_over_dir,'AP')
        FOV = [DWIparams.FOV.FH, DWIparams.FOV.AP]; % [M,P]
    elseif strcmpi(DWIparams.fold_over_dir,'FH')
        FOV = [DWIparams.FOV.AP, DWIparams.FOV.FH]; % [M,P]
    else
        error('f_get_recon_params:main','Unknown DWI fold-over direction')
    end
elseif strcmpi(DWIparams.slice_orientation,'coronal')    
    if strcmpi(DWIparams.fold_over_dir,'RL')
        FOV = [DWIparams.FOV.FH, DWIparams.FOV.RL]; % [M,P]
    elseif strcmpi(DWIparams.fold_over_dir,'FH')
        FOV = [DWIparams.FOV.RL, DWIparams.FOV.FH]; % [M,P]
    else
        error('f_get_recon_params:main','Unknown DWI fold-over direction')
    end
else
    error('f_get_recon_params:main','Unknown DWI slice orientation')
end
UGN1_ACQ_reduced_acq_time = DWIparams.UGN1_ACQ_reduced_acq_time;
% rec_vox_p = DWIparams.acq_vox_m; % recon isotropic in-plane voxel size

M_ORI = 1; % readout (measurement) direction
P_ORI = 2; % phase-encoding direction



%% Calculate recon matrix and voxel size for SMAP

% Get scan parameters. 
%* Part of these can be verified from ExamCard.
cd(DWIparams.loadDir_s) % to read .LIST file in [f_get_params_from_listfile.m]
[ENCima,LCAima,STACKima,ACQREC] = f_get_ENCima_v2(ACQ_RESOL,RECON_RESOL, ...
    SENSE_FACTOR,isEPI,EPI_FACTOR1,slice_thickness,FOV, ...
    UGN1_ACQ_reduced_acq_time,STACKima.aspect_ratio,DWIparams);


% Find final image matrix size of SMAP along the phase-encoding and
% frequency-encoding direction of ALIASED data.

% Get the size for resizing the SMAP matrix along P_ORI and M_ORI of
% ALIASED data. This is to match the SMAP voxel size to ALIASED data voxel
% size along the P_ORI and M_ORI of ALIASED data.

% Assumption is that the REF and DWI is acquired on the same slice
% orientation.

if strcmpi(DWIparams.slice_orientation,'transverse')
    
    if strcmpi(DWIparams.fold_over_dir,'RL')                
        if strcmpi(REFparams.fold_over_dir,'AP')            
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP M_ORI.
            s_FOV_p = REFparams.FOV.RL; % SMAP
            s_acq_vox_p = REFparams.acq_vox.RL;
            s_ovs_fac_p = REFparams.ovs_fac_m;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP P_ORI.
            s_FOV_m = REFparams.FOV.AP; % SMAP
            s_acq_vox_m = REFparams.acq_vox.AP;
            s_ovs_fac_m = REFparams.ovs_fac_p;            
        elseif strcmpi(REFparams.fold_over_dir,'RL')            
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP P_ORI.
            s_FOV_p = REFparams.FOV.RL; % SMAP
            s_acq_vox_p = REFparams.acq_vox.RL;
            s_ovs_fac_p = REFparams.ovs_fac_p;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP M_ORI.
            s_FOV_m = REFparams.FOV.AP; % SMAP
            s_acq_vox_m = REFparams.acq_vox.AP;
            s_ovs_fac_m = REFparams.ovs_fac_m;
        else
            error('f_get_recon_params:main','Unknown REFparams.fold_over_dir')
        end
        %---------- P_ORI ----------        
        a_FOV_p = DWIparams.FOV.RL; % ALIASED (image-echo)
        a_acq_vox_p = ACQREC.ACQ_voxel_MPS(P_ORI);
        a_ovs_fac_p = ENCima.oversample_factors(P_ORI);                
        %---------- M_ORI ----------
        a_FOV_m = DWIparams.FOV.AP; % ALIASED (image-echo)
        a_acq_vox_m = ACQREC.ACQ_voxel_MPS(M_ORI);
        a_ovs_fac_m = ENCima.oversample_factors(M_ORI);
        
    elseif strcmpi(DWIparams.fold_over_dir,'AP')        
        if strcmpi(REFparams.fold_over_dir,'AP')
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP P_ORI.
            s_FOV_p = REFparams.FOV.AP;
            s_acq_vox_p = REFparams.acq_vox.AP;
            s_ovs_fac_p = REFparams.ovs_fac_p;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP M_ORI.
            s_FOV_m = REFparams.FOV.RL; % SMAP
            s_acq_vox_m = REFparams.acq_vox.RL;
            s_ovs_fac_m = REFparams.ovs_fac_m; 
            
        elseif strcmpi(REFparams.fold_over_dir,'RL')
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP M_ORI.
            s_FOV_p = REFparams.FOV.AP;
            s_acq_vox_p = REFparams.acq_vox.AP;
            s_ovs_fac_p = REFparams.ovs_fac_m;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP P_ORI.
            s_FOV_m = REFparams.FOV.RL; % SMAP
            s_acq_vox_m = REFparams.acq_vox.RL;
            s_ovs_fac_m = REFparams.ovs_fac_p;
            
        elseif isempty(REFparams.fold_over_dir)
            s_FOV_p = [];
            s_acq_vox_p = [];
            s_ovs_fac_p = [];
            s_FOV_m = [];
            s_acq_vox_m = [];
            s_ovs_fac_m = [];
        else
            error('f_get_recon_params:main','Unknown REFparams.fold_over_dir')
        end
        %---------- P_ORI ----------
        a_FOV_p = DWIparams.FOV.AP;
        a_acq_vox_p = ACQREC.ACQ_voxel_MPS(P_ORI);
        a_ovs_fac_p = ENCima.oversample_factors(P_ORI);
        %---------- M_ORI ----------
        a_FOV_m = DWIparams.FOV.RL; % ALIASED (image-echo)
        a_acq_vox_m = ACQREC.ACQ_voxel_MPS(M_ORI);
        a_ovs_fac_m = ENCima.oversample_factors(M_ORI);
        
    else
        error('f_get_recon_params:main','Unknown DWI fold-over direction')
    end
    
elseif strcmpi(DWIparams.slice_orientation,'sagittal')
    
    if strcmpi(DWIparams.fold_over_dir,'AP')                
        if strcmpi(REFparams.fold_over_dir,'AP')
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP P_ORI.
            s_FOV_p = REFparams.FOV.AP; % SMAP
            s_acq_vox_p = REFparams.acq_vox.AP;
            s_ovs_fac_p = REFparams.ovs_fac_p;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP M_ORI.
            s_FOV_m = REFparams.FOV.FH; % SMAP
            s_acq_vox_m = REFparams.acq_vox.FH;
            s_ovs_fac_m = REFparams.ovs_fac_m;
            
        elseif strcmpi(REFparams.fold_over_dir,'FH')
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP M_ORI.
            s_FOV_p = REFparams.FOV.AP; % SMAP
            s_acq_vox_p = REFparams.acq_vox.AP;
            s_ovs_fac_p = REFparams.ovs_fac_m;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP P_ORI.
            s_FOV_m = REFparams.FOV.FH; % SMAP
            s_acq_vox_m = REFparams.acq_vox.FH;
            s_ovs_fac_m = REFparams.ovs_fac_p;
            
        elseif isempty(REFparams.fold_over_dir)
            s_FOV_p = [];
            s_acq_vox_p = [];
            s_ovs_fac_p = [];
            s_FOV_m = [];
            s_acq_vox_m = [];
            s_ovs_fac_m = [];
            
        else
            error('f_get_recon_params:main','Unknown REFparams.fold_over_dir')
        end
        %---------- P_ORI ----------
        a_FOV_p = DWIparams.FOV.AP; % ALIASED (image-echo)
        a_acq_vox_p = ACQREC.ACQ_voxel_MPS(P_ORI);
        a_ovs_fac_p = ENCima.oversample_factors(P_ORI);
        %---------- M_ORI ----------
        a_FOV_m = DWIparams.FOV.FH; % ALIASED (image-echo)
        a_acq_vox_m = ACQREC.ACQ_voxel_MPS(M_ORI);
        a_ovs_fac_m = ENCima.oversample_factors(M_ORI);
        
    elseif strcmpi(DWIparams.fold_over_dir,'FH')
        if strcmpi(REFparams.fold_over_dir,'FH')
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP P_ORI.
            s_FOV_p = REFparams.FOV.FH;
            s_acq_vox_p = REFparams.acq_vox.FH;
            s_ovs_fac_p = REFparams.ovs_fac_p;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP M_ORI.
            s_FOV_m = REFparams.FOV.AP; % SMAP
            s_acq_vox_m = REFparams.acq_vox.AP;
            s_ovs_fac_m = REFparams.ovs_fac_m;
            
        elseif strcmpi(REFparams.fold_over_dir,'AP')
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP M_ORI.
            s_FOV_p = REFparams.FOV.FH;
            s_acq_vox_p = REFparams.acq_vox.FH;
            s_ovs_fac_p = REFparams.ovs_fac_m;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP P_ORI.
            s_FOV_m = REFparams.FOV.AP; % SMAP
            s_acq_vox_m = REFparams.acq_vox.AP;
            s_ovs_fac_m = REFparams.ovs_fac_p;
            
        elseif isempty(REFparams.fold_over_dir)
            s_FOV_p = [];
            s_acq_vox_p = [];
            s_ovs_fac_p = [];
            s_FOV_m = [];
            s_acq_vox_m = [];
            s_ovs_fac_m = [];
            
        else
            error('f_get_recon_params:main','Unknown REFparams.fold_over_dir')
        end
        %---------- P_ORI ----------
        a_FOV_p = DWIparams.FOV.FH;
        a_acq_vox_p = ACQREC.ACQ_voxel_MPS(P_ORI);
        a_ovs_fac_p = ENCima.oversample_factors(P_ORI);
        %---------- M_ORI ----------
        a_FOV_m = DWIparams.FOV.AP; % ALIASED (image-echo)
        a_acq_vox_m = ACQREC.ACQ_voxel_MPS(M_ORI);
        a_ovs_fac_m = ENCima.oversample_factors(M_ORI);
        
    else
        error('f_get_recon_params:main','Unknown DWI fold-over direction')
    end
    
elseif strcmpi(DWIparams.slice_orientation,'coronal')
    
    if strcmpi(DWIparams.fold_over_dir,'FH')            
        if strcmpi(REFparams.fold_over_dir,'FH')
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP P_ORI.
            s_FOV_p = REFparams.FOV.FH; % SMAP
            s_acq_vox_p = REFparams.acq_vox.FH;
            s_ovs_fac_p = REFparams.ovs_fac_p;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP M_ORI.
            s_FOV_m = REFparams.FOV.RL; % SMAP
            s_acq_vox_m = REFparams.acq_vox.RL;
            s_ovs_fac_m = REFparams.ovs_fac_m;
            
        elseif strcmpi(REFparams.fold_over_dir,'RL')
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP M_ORI.
            s_FOV_p = REFparams.FOV.FH; % SMAP
            s_acq_vox_p = REFparams.acq_vox.FH;
            s_ovs_fac_p = REFparams.ovs_fac_m;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP P_ORI.
            s_FOV_m = REFparams.FOV.RL; % SMAP
            s_acq_vox_m = REFparams.acq_vox.RL;
            s_ovs_fac_m = REFparams.ovs_fac_p;
            
        elseif isempty(REFparams.fold_over_dir)
            s_FOV_p = [];
            s_acq_vox_p = [];
            s_ovs_fac_p = [];
            s_FOV_m = [];
            s_acq_vox_m = [];
            s_ovs_fac_m = [];
            
        else
            error('f_get_recon_params:main','Unknown REFparams.fold_over_dir')
        end
        %---------- P_ORI ----------
        a_FOV_p = DWIparams.FOV.FH; % ALIASED (image-echo)
        a_acq_vox_p = ACQREC.ACQ_voxel_MPS(P_ORI);
        a_ovs_fac_p = ENCima.oversample_factors(P_ORI);
        %---------- M_ORI ----------
        a_FOV_m = DWIparams.FOV.RL; % ALIASED (image-echo)
        a_acq_vox_m = ACQREC.ACQ_voxel_MPS(M_ORI);
        a_ovs_fac_m = ENCima.oversample_factors(M_ORI);
        
    elseif strcmpi(DWIparams.fold_over_dir,'RL')
        if strcmpi(REFparams.fold_over_dir,'RL')
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP P_ORI.
            s_FOV_p = REFparams.FOV.RL;
            s_acq_vox_p = REFparams.acq_vox.RL;
            s_ovs_fac_p = REFparams.ovs_fac_p;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP M_ORI.
            s_FOV_m = REFparams.FOV.FH; % SMAP
            s_acq_vox_m = REFparams.acq_vox.FH;
            s_ovs_fac_m = REFparams.ovs_fac_m;
            
        elseif strcmpi(REFparams.fold_over_dir,'FH')
            %---------- P_ORI ----------
            % SMAP in ALIASED P_ORI == SMAP M_ORI.
            s_FOV_p = REFparams.FOV.RL;
            s_acq_vox_p = REFparams.acq_vox.RL;
            s_ovs_fac_p = REFparams.ovs_fac_m;
            
            %---------- M_ORI ----------
            % SMAP ALIASED M_ORI == SMAP P_ORI.
            s_FOV_m = REFparams.FOV.FH; % SMAP
            s_acq_vox_m = REFparams.acq_vox.FH;
            s_ovs_fac_m = REFparams.ovs_fac_p;
            
        elseif isempty(REFparams.fold_over_dir)
            s_FOV_p = [];
            s_acq_vox_p = [];
            s_ovs_fac_p = [];
            s_FOV_m = [];
            s_acq_vox_m = [];
            s_ovs_fac_m = [];
            
        else
            error('f_get_recon_params:main','Unknown REFparams.fold_over_dir')
        end
        %---------- P_ORI ----------
        a_FOV_p = DWIparams.FOV.RL;
        a_acq_vox_p = ACQREC.ACQ_voxel_MPS(P_ORI);
        a_ovs_fac_p = ENCima.oversample_factors(P_ORI);
        %---------- M_ORI ----------
        a_FOV_m = DWIparams.FOV.FH; % ALIASED (image-echo)
        a_acq_vox_m = ACQREC.ACQ_voxel_MPS(M_ORI);
        a_ovs_fac_m = ENCima.oversample_factors(M_ORI);
        
    else
        error('f_get_recon_params:main','Unknown DWI fold-over direction')
    end
else
    error('f_get_recon_params:main','Unknown DWI slice orientation')
end


%-------------------- P_ORI --------------------
% Matrix size N of SMAP along P_ORI of ALIASED data.
s_N_p = s_FOV_p / s_acq_vox_p * s_ovs_fac_p;


% Actual matrix size N of ALIASED after SENSE recon along P_ORI.
%* 'ACQ matrix M x P' in ExamCard is not actual acquisition matrix for EPI case.
%* This is (epi_fac * shots * sense_fac).
%* (FOV / acq_vox_p * ovs_fac_p = epi_fac * shots * sense_fac = ovs_res * sense_fac).
a_N_p = a_FOV_p /  a_acq_vox_p * a_ovs_fac_p;


% Final recon size of ALIASED along P_ORI.
E = DWIparams.EPI_FACTOR(1);
ns = ACQREC.number_of_shots;
R = double(uint8(ENCima.sense_factors(P_ORI))); % 2010.08.17
% a_N_recon_p = E * ns * R; % a_N_p must be the same as a_N_recon_p

% Above can generate problem in '% Determine if s_Np_res must be zeropadded or cut.' part below.
%* 2010.07.22.
a_N_recon_p = E * ns * R; 


% Get the SMAP matrix size to be resized along P_ORI of ALIASED data.
%* After resize along P_ORI of ALIASED data, SMAP voxel size along P_ORI of 
%* ALIASED data is the same as the voxel size of ALIASED data along P_ORI.
s_N_res_p = s_FOV_p / a_acq_vox_p * s_ovs_fac_p;


% Determine if s_Np_res must be zeropadded or cut.
s_cut_vox_p = [];
s_zeropad_vox_p = [];
if s_N_res_p > a_N_recon_p % cut
    if mod(ceil(s_N_res_p) - a_N_recon_p,2)==0
        s_cut_vox_p = ( ceil(s_N_res_p) - a_N_recon_p )/2;
        s_N_res_p = ceil(s_N_res_p);
    else
        s_cut_vox_p = ( floor(s_N_res_p) - a_N_recon_p )/2;
        s_N_res_p = floor(s_N_res_p);
    end
else % zeropad
    if mod(a_N_recon_p - ceil(s_N_res_p),2)==0
        s_zeropad_vox_p = ( a_N_recon_p - ceil(s_N_res_p) )/2;
        s_N_res_p = ceil(s_N_res_p);
    else
        s_zeropad_vox_p = ( a_N_recon_p - floor(s_N_res_p) )/2;
        s_N_res_p = floor(s_N_res_p);
    end
end


% Final matrix size of SMAP along P_ORI of ALIASED data.
s_N_final_p = a_N_recon_p; % same as a_N_p


% CHECK IF,
%a_N_p == a_N_recon_p == E * ns * R == s_N_final_p



%-------------------- M_ORI --------------------
% Matrix size N of SMAP along M_ORI of ALIASED data.
s_N_m = s_FOV_m / s_acq_vox_m * s_ovs_fac_m;


% Actual matrix size N of ALIASED after SENSE recon along M_ORI.
%* 'ACQ matrix M x P' in ExamCard is not actual acquisition matrix for EPI case.
%* This is (epi_fac * shots * sense_fac).
%* (FOV / acq_vox_m * ovs_fac_m = epi_fac * shots * sense_fac = ovs_res * sense_fac).
a_N_m = a_FOV_m /  a_acq_vox_m * a_ovs_fac_m;


% Final recon size of ALIASED along M_ORI.
a_N_recon_m = a_N_m; % a_N_m must be the same as a_N_recon_m


% Get the SMAP matrix size to be resized along M_ORI of ALIASED data.
%* After resize along M_ORI of ALIASED data, SMAP voxel size along M_ORI of 
%* ALIASED data is the same as the voxel size of ALIASED data along M_ORI.
s_N_res_m = s_FOV_m / a_acq_vox_m * s_ovs_fac_m;


% Determine if s_Np_res must be zeropadded or cut.
s_cut_vox_m = [];
s_zeropad_vox_m = [];
if s_N_res_m > a_N_recon_m % cut
    if mod(ceil(s_N_res_m) - a_N_recon_m,2)==0
        s_cut_vox_m = ( ceil(s_N_res_m) - a_N_recon_m )/2;
        s_N_res_m = ceil(s_N_res_m);
    else
        s_cut_vox_m = ( floor(s_N_res_m) - a_N_recon_m )/2;
        s_N_res_m = floor(s_N_res_m);
    end
else % zeropad
    if mod(a_N_recon_m - ceil(s_N_res_m),2)==0
        s_zeropad_vox_m = ( a_N_recon_m - ceil(s_N_res_m) )/2;
        s_N_res_m = ceil(s_N_res_m);
    else
        s_zeropad_vox_m = ( a_N_recon_m - floor(s_N_res_m) )/2;
        s_N_res_m = floor(s_N_res_m);
    end
end


% Final matrix size of SMAP along M_ORI of ALIASED data.
s_N_final_m = a_N_recon_m; % same as a_N_m


% CHECK IF,
%a_N_m == a_N_recon_m == s_N_fianl_m



%% Calculate final recon matrix for ALIASED

% Acquisition voxel size is
% acq_vox_p = FOV * ovs_fac_p / (epi_fac * ns * sense_fac)
%
% Then the number of voxels to have recon_vox_p size voxel along P_ORI is
% acq_vox_p/ovs_fac_p : FOV/(epi_fac * ns * sense_fac) == recon_vox_p : FOV / N_recon
%
% here,
% recon_vox_p: desired recon voxel size for FOV, which is using interpolation after recon
% N_recon: number of voxels corresponding to recon_vox_p for the same FOV
%
% then,
% (acq_vox_p/ovs_fac_p) : 1/(epi_fac*ns*sense_fac) = recon_vox_p : (1/N_recon).
%
% N_recon = acq_vox_p * (epi_fac*ns*sense_fac) / (ovs_fac_p * recon_vox_p).


% This is the matrix size to have ACQREC.REC_voxel_MPS(M_ORI) voxel size 
% along P_ORI.

a_N_final_interp_p = ACQREC.ACQ_voxel_MPS(P_ORI) * ...
    (ENCima.epi_factors * ACQREC.number_of_shots * ENCima.sense_factors(P_ORI)) / ...
    (ENCima.oversample_factors(P_ORI) * ACQREC.REC_voxel_MPS(M_ORI)); % same as DWIparams.FOV(P_ORI)
% a_N_final_interp_p = ACQREC.ACQ_voxel_MPS(P_ORI) * ...
%     (ENCima.epi_factors * ACQREC.number_of_shots * ENCima.sense_factors(P_ORI)) / ...
%     (ACQREC.REC_voxel_MPS(M_ORI));

%* The reason why this use ACQREC.REC_voxel_MPS(M_ORI) is that we try to recon
%* in isotropic in-plane voxels.


% This is the matrix size after zeropad to fill to have the square image matrix.
a_N_final_interp_recon_p = ENCima.recon_resolutions(M_ORI);
if mod(a_N_final_interp_recon_p - ceil(a_N_final_interp_p),2)==0
    a_zeropad_vox_p = (a_N_final_interp_recon_p - ceil(a_N_final_interp_p))/2;
else
    a_zeropad_vox_p = (a_N_final_interp_recon_p - floor(a_N_final_interp_p))/2;
end



%% Set RECONparams
clear  RECONparams

RECONparams.ENCima = ENCima;
RECONparams.LCAima = LCAima;
RECONparams.STACKima = STACKima;
RECONparams.ACQREC = ACQREC;

RECONparams.s_N_res_m = s_N_res_m;
RECONparams.s_N_final_m = s_N_final_m;
RECONparams.s_N_res_p = s_N_res_p;
RECONparams.s_N_final_p = s_N_final_p;

% RECONparams.a_N_recon_m = a_N_recon_m;
% RECONparams.a_N_recon_p = a_N_recon_p;

RECONparams.s_cut_vox_m = s_cut_vox_m;
RECONparams.s_zeropad_vox_m = s_zeropad_vox_m;
RECONparams.s_cut_vox_p = s_cut_vox_p;
RECONparams.s_zeropad_vox_p = s_zeropad_vox_p;

RECONparams.a_N_recon_m = a_N_recon_m;
RECONparams.a_N_recon_p = a_N_recon_p;
RECONparams.a_N_final_interp_p = a_N_final_interp_p;
RECONparams.a_N_final_interp_recon_p = a_N_final_interp_recon_p;
RECONparams.a_zeropad_vox_p = a_zeropad_vox_p;



%% Report
%----- This is taken care in [s_check_smapSizeInfo.m] -----
% if verbose==1
%     fprintf('    SMAP adjusting parameters for IMG recon.\n')
%     fprintf('       SMAP    FOV                     : [M,P] = [%d,%d]\n', ...
%         REFparams.FOV(2),REFparams.FOV(1))
%     fprintf('       SMAP    acq_vox                 : [M,P] = [%.3f,%.3f]\n', ...
%         REFparams.acq_vox_m,REFparams.acq_vox_p)
%     fprintf('       SMAP    ovs_fac                 : [M,P] = [%.3f,%.3f]\n', ...
%         REFparams.ovs_fac_m,REFparams.ovs_fac_p)
%     fprintf('       SMAP    matrix                  : [M,P] = [%.3f,%.3f]\n', ...
%         [REFparams.FOV(2),REFparams.FOV(1)] ./ [REFparams.acq_vox_m,REFparams.acq_vox_p] .* ...
%         [REFparams.ovs_fac_m,REFparams.ovs_fac_p])
%     
%     disp(' ')
%     
%     fprintf('       ALIASED FOV                     : [M,P] = [%d,%d]\n', ...
%         DWIparams.FOV(1),DWIparams.FOV(2))
%     fprintf('       ALIASED acq_vox                 : [M,P] = [%.3f,%.3f]\n', ...
%         RECONparams.ACQREC.ACQ_voxel_MPS(1:2))
%     fprintf('       ALIASED ovs_fac                 : [M,P] = [%.3f,%.3f]\n', ...
%         RECONparams.ENCima.oversample_factors(1:2))
%     fprintf('       ALIASED matrix (after recon)    : [M,P] = [%.3f,%.3f]\n', ...
%         [DWIparams.FOV(1),DWIparams.FOV(2)] ./ RECONparams.ACQREC.ACQ_voxel_MPS(1:2) .* ...
%         RECONparams.ENCima.oversample_factors(1:2))
%     
%     disp(' ')
%     
%     fprintf('       SMAP    resize                  : [M,P] = [%.3f,%.3f]\n', ...
%         RECONparams.s_N_res_m, RECONparams.s_N_res_p)
%     fprintf('       SMAP    cut                     : [M,P] = [%.3f,%.3f]\n', ...
%         RECONparams.s_cut_vox_m, RECONparams.s_cut_vox_p)
%     fprintf('       SMAP    zeropad                 : [M,P] = [%.3f,%.3f]\n', ...
%         RECONparams.s_zeropad_vox_m, RECONparams.s_zeropad_vox_p)
%     fprintf('       SMAP    after cut and zeropad   : [M,P] = [%.3f,%.3f]\n', ...
%         RECONparams.s_N_final_m, RECONparams.s_N_final_p)
%     
%     fprintf('\n')
% end



%% END










