
%[s_get_REFparams] get coil reference scan parameters, REFparams.
%
% USAGE:
%   s_get_REFparams
%
% OUTPUT:
%   REFparams: A struct with fields as,
%       filename
%       FOV, image field-of-view, [AP,RL]
%       acq_vox_m, acquisition voxel size along readout
%       acq_vox_p, acquisition voxel size along phase-encoding
%       slice_thickness
%       ovs_fac_m, oversampling factor along readout 
%       ovs_fac_p, oversampling factor along phase-encoding 
%       
%
%
% Last modified
% 2010.08.06.
% 2010.09.01.
%   REF FOV has fields [AP,RL,FH].
%   REF acq_vox has fields [AP,RL,FH].
% 2010.09.09.
%   Add REFparams.bodycoil_idx. This is bodycoil index in the .LIST file
%   starting from 0.
% 2010.09.16.
%   Add 'SENSE-Breast-4' coil and REFparams.slice_orientation.
% 2010.09.23.
%   Add fat_shift_dir.
% 2010.09.24.
%   Adjust ovs_fac_m(p) for 'SENSE-NV-16' coil.
% 2010.10.08.
%   Adjust ovs_fac_p when 'Fold-over suppression' = 'yes' from 1 to 2.
% 2010.12.06.
%   Adjust REFparams.slice_thickness based on refparams.slice_orientation.
% 2011.03.03.
%   Add Offc. parameters.
% 2011.03.31.
%   Add SENSE-Torso coil for 3T.
% 2012.02.24.
%   Add SENSE-Head-32 coil for 3T.
% 2012.04.16.
%   Add RX-Intf-1_Quad-TR-1 for 7T 16 ch Spine coil.
% 2012.05.22.
%   Take care of 'no SENSE' case.
%
% Ha-Kyu



%% Define REFparams

% Clear.
clear  REFparams


% Filename.
REFparams.filename = REFrawfilename;

if ~isempty(REFparams.filename)
    
% GEOMETRY.
REFparams.FOV = refparams.fov; % refparams.fov.[AP,RL,FH]
REFparams.acq_vox = refparams.acq_vox; % refparams.acq_vox.[AP,RL,FH]
if strcmpi(refparams.slice_orientation,'transverse')
    REFparams.slice_thickness = refparams.acq_vox.FH;
elseif strcmpi(refparams.slice_orientation,'coronal')
    REFparams.slice_thickness = refparams.acq_vox.AP;
elseif strcmpi(refparams.slice_orientation,'sagittal')
    REFparams.slice_thickness = refparams.acq_vox.RL;
else
    error('s_get_REFparams:main','Unknown refparams.slice_orientation')
end
REFparams.fold_over_suppression = refparams.fold_over_suppression; % yes or no
REFparams.slices = refparams.slices;
REFparams.slice_thickness = refparams.slice_thickness;
REFparams.slice_orientation = refparams.slice_orientation; 
REFparams.fold_over_dir = refparams.fold_over_dir;
REFparams.fat_shift_dir = refparams.fat_shift_dir;


% OFFC/ANG.
REFparams.Offc.AP = refparams.offc(1);
REFparams.Offc.RL = refparams.offc(2);
REFparams.Offc.FH = refparams.offc(3);
REFparams.Ang.AP = refparams.ang(1);
REFparams.Ang.RL = refparams.ang(2);
REFparams.Ang.FH = refparams.ang(3);


% Internal.
REFparams.ovs_fac_m = 2.0; % default:2 for readout
REFparams.ovs_fac_p = 1.0; % default:1 for gradient echo
if strcmpi(REFparams.fold_over_suppression,'yes')
    REFparams.ovs_fac_p = 2.0; % 2.0 when 'Fold-over suppression = yes'
end


% Bodycoil index (starts from 0).
if GENparams.B0==70000    
    if strcmpi(GENparams.coilID,'SENSE-Head-7T')
        REFparams.nCOIL = 16;
        REFparams.bodycoil_idx = 26; % starting from 0
    elseif strcmpi(GENparams.coilID,'RX-Intf-1')
        REFparams.nCOIL = 32;
        REFparams.bodycoil_idx = 26; % starting from 0
    elseif strcmpi(GENparams.coilID,'RX-Intf-1_Quad-TR-1')
        REFparams.nCOIL = 16;
        REFparams.bodycoil_idx = 26; % starting from 0
    else
        error('s_get_REFparams:main','Other GENparams.coilID is not yet defined')
    end
elseif GENparams.B0==30000    
    if strcmpi(GENparams.coilID,'SENSE-Head-8')
        REFparams.nCOIL = 8;
        REFparams.bodycoil_idx = 27; % starting from 0
    elseif strcmpi(GENparams.coilID,'SENSE-NV-16')
        REFparams.nCOIL = 16;
        %REFparams.nCOIL = 6; % H elemens selection
        REFparams.bodycoil_idx = 27; % starting from 0 --> need to check as of 2010.09.09
    elseif strcmpi(GENparams.coilID,'SENSE-Breast-4')
        REFparams.nCOIL = 4;
        REFparams.bodycoil_idx = 27;
    elseif strcmpi(GENparams.coilID,'SENSE-Torso')
        REFparams.nCOIL = 6;
        REFparams.bodycoil_idx = 27;
    elseif strcmpi(GENparams.coilID,'SENSE-Head-32')
        REFparams.nCOIL = 32;
        REFparams.bodycoil_idx = 27;
    else
        error('s_get_REFparams:main','Other GENparams.coilID is not yet defined')
    end
else
    error('s_get_REFparams:main','Unknown GENparams.B0')
end

else
    
    REFparams.FOV = [];
    REFparams.acq_vox = [];
    REFparams.slice_thickness = [];
    REFparams.fold_over_suppression = [];
    REFparams.slice_orientation = [];
    REFparams.fold_over_dir = [];
    REFparams.fat_shift_dir = [];
    REFparams.Offc.AP = [];
    REFparams.Offc.RL = [];
    REFparams.Offc.FH = [];
    REFparams.Ang.AP = [];
    REFparams.Ang.RL = [];
    REFparams.Ang.FH = [];
    REFparams.ovs_fac_m = [];
    REFparams.ovs_fac_p = [];
    REFparams.nCOIL = [];
    REFparams.bodycoil_idx = [];
        
end





