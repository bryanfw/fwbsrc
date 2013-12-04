
%[s_get_FMAPparams] get FMAP parameters, FMAPparams.
%
% USAGE:
%   s_get_FMAPparams
%
% OUTPUT:
%   FMAPparams: A struct with fields as,
%       acq_vox_m, ACQ voxel MPS (M_ORI)
%       acq_vox_p, ACQ voxel MPS (P_ORI)
%       slice_thickness, ACQ voxel MPS (S_ORI)
%       rec_vox_m, Reconstruction voxel size (M_ORI)
%       rec_vox_p, Reconstruction voxel size (P_ORI)
%       filename, Fieldmap PAR/REF filename
%       fitorder, Order of 2D polynomial fitting
%       fitsize, [fitsize_out,fitsize_in], local 2D polynomial fitting of fieldmap
%       vox_dilate, Number of voxels to dilate B0 fieldmap mask
%       vox_kern, Kernel size (voxels) for smoothing B0 fieldmap
%       fitregion, [local,global,none] for local, global and no 2D polynomial
%
%
% Last modified
% 2010.08.06.
% 2010.08.13.
%   Use [f_read_parrec.m].
% 2010.09.16.
%   Modify for sagittal slice orientation.
% 2010.09.22.
%   Get FMAPparams.filename even in empty filename case.
% 2011.03.28.
%   Modify FMAPparams.vox_kern.



%% Define FMAPparams

if ~isempty(FMAPparrecfilename)
    clear  FMAPparams
    
    [data,info] = f_read_parrec(fullfile(loadParDir_s,FMAPparrecfilename));
    FMAPparams.filename =  FMAPparrecfilename;
    scan_res_v = str2num(info.Scan_resolution); % [M,P]
    fov_v = str2num(info.FOV); % [ap,fh,rl]
    
    % Reserve FOV.
    FMAPparams.FOV.AP = fov_v(1);
    FMAPparams.FOV.FH = fov_v(2);
    FMAPparams.FOV.RL = fov_v(3);
    
    % Take FOV [M,P] and acq_vox [M,P].
    fold_over_dir = info.Preparation_direction;
    if (info.imgdef.slice_orientation.uniq==1) %[1,2,3]==[transverse,sagittal,coronal]
        if strcmpi(fold_over_dir,'Anterior-Posterior')
            fov_v = fov_v([3,1]); % must be [M,P]==[rl,ap]
            fold_over_dir = 'AP';
        elseif strcmpi(fold_over_dir,'Right-left')
            fov_v = fov_v([1,3]); % [M,P]==[ap,rl]
            fold_over_dir = 'RL';
        else
            error('s_get_FMAPparams:main','Unsupported FMAP info.Preparation_direction')
        end
        slice_orientation = 'transverse';
    elseif (info.imgdef.slice_orientation.uniq==2)
        if strcmpi(fold_over_dir,'Anterior-Posterior')
            fov_v = fov_v([2,1]); % [M,P]==[fh,ap]
            fold_over_dir = 'AP';
        elseif strcmpi(fold_over_dir,'Feet-Head') %is Feet-Head correct?
            fov_v = fov_v([1,2]); % [M,P]==[ap,fh]
            fold_over_dir = 'FH';
        else
            error('s_get_FMAPparams:main','Unsupported FMAP info.Preparation_direction')
        end
        slice_orientation = 'sagittal';
    elseif (info.imgdef.slice_orientation.uniq==3)
        if strcmpi(fold_over_dir,'Right-left')
            fov_v = fov_v([2,3]); % [M,P]==[fh,rl]
            fold_over_dir = 'RL';
        elseif strcmpi(fold_over_dir,'Feet-Head')
            fov_v = fov_v([3,2]); % [M,P]==[rl,fh]
            fold_over_dir = 'FH';
        else
            error('s_get_FMAPparams:main','Unsupported FMAP info.Preparation_direction')
        end
        slice_orientation = 'coronal';
    else
        error('s_get_FMAPparams:main','Unknown FMAP info.imgdef.slice_orientation.uniq')
    end
    acq_vox_v = fov_v./scan_res_v;
    %if strcmpi(info.Preparation_direction,'Right-left')
    %    FMAPparams.acq_vox_m = acq_vox_v(1); % ACQ voxel MPS (M_ORI)
    %    FMAPparams.acq_vox_p = acq_vox_v(2); % ACQ voxel MPS (P_ORI)
    %else
    %    FMAPparams.acq_vox_m = acq_vox_v(2);
    %    FMAPparams.acq_vox_p = acq_vox_v(1);
    %end
    
    % Reserve other parameters.
    FMAPparams.acq_vox_m = acq_vox_v(1); % always M_ORI
    FMAPparams.acq_vox_p = acq_vox_v(2); % always P_ORI
    FMAPparams.slice_thickness = info.imgdef.slice_thickness.uniq; % ACQ voxel MPS (S_ORI)
    rec_vox_mp = info.imgdef.pixel_spacing.uniq;
    FMAPparams.rec_vox_m = rec_vox_mp(1); % Reconstruction voxel size (M_ORI)
    FMAPparams.rec_vox_p = rec_vox_mp(2); % Reconstruction voxel size (P_ORI)
    FMAPparams.slice_orientation = slice_orientation; % [transverse,sagittal,coronal]
    FMAPparams.fold_over_dir = fold_over_dir; % fold-over direction
else
    FMAPparams.filename =  FMAPparrecfilename; % empty
    FMAPparams.acq_vox_m = [];
    FMAPparams.acq_vox_p = [];
    FMAPparams.slice_thickness = [];
    FMAPparams.rec_vox_m = [];
    FMAPparams.rec_vox_m = [];
    FMAPparams.slice_orientation = [];
    FMAPparams.fold_over_dir = [];
end


% Fit paramters.
fitorder = 8;       % order of 2D polynomial fitting, (default:8)
% vox_kern = 3;       % kernel size (voxels) for smoothing B0 fieldmap, (default:3)
vox_dilate = 3;     % number of voxels to dilate B0 fieldmap mask, (default:0)
fitsize_out = 5;    % for local 2D polynomial fitting of fieldmap (outside object), (default:5)
fitsize_in = 3;     % for local 2D polynomial fitting of fieldmap (inside object), (default:3)

FMAPparams.fitorder = fitorder;
FMAPparams.fitsize = [fitsize_out,fitsize_in];
FMAPparams.vox_dilate = vox_dilate;
% FMAPparams.vox_kern = RECONparams.ACQREC.ACQ_voxel_MPS(2) * RECONparams.ACQREC.number_of_shots * 2;
FMAPparams.vox_kern = RECONparams.ACQREC.ACQ_voxel_MPS(2) * 2;

% Fit region; [local,global,none] for local, global and no 2D polynomial
% fitting
FMAPparams.fitregion = 'none';%'global'; 

clear  fitorder  fitsize_out  fitsize_in  vox_dilate  vox_kern  slice_orientation











