
function [CRP,NWV,MPS,LPH,LPHp,XYZ] = f_coord_tform(PHL, Ang, ...
    fold_over_dir, fat_shift_dir, slice_orientation, ...
    patient_orientation, patient_position)
%[f_coord_tform] transform scanner coordinate from PHL to all the other.
%
% USAGE:
% [CRP,NWV,MPS,LPH,LPHp,XYZ] = f_coord_tform(PHL, Ang, ...
%     fold_over_dir, fat_shift_dir, slice_orientation, ...
%     patient_orientation, patient_position)
%
% INPUT:
%   PHL - 
%       Scanner coordinate vector in P(AP)-H(FH)-L(RL) order in PHL
%       coordinate system. This is the same order as shown in .PAR file for
%       diffusion direction.
%   Ang - 
%       Slice angulation with fields,
%           AP, in degrees
%           RL, in degrees
%           FL, in degrees
%   fold_over_dir -
%       Fold-over direction (phase-encoding direction)
%   fat_shift_dir -
%       Fat-shift direction
%   slice_orientation -
%       Slice orientation (transverse, coronal, sagittal)
%   patient_orientation - 
%       supine, prone, rightdecubitus or leftdecubitus
%   patient_position -
%       headfirst or feetfirst
%
% OUTPUT:
%   CRP - vectors in CRP coordinate
%   NWV - vectors in NWV coordinate
%   MPS - vectors in MPS coordinate
%   LPH - vectors in LPH coordinate
%   LPHp - vectors in LPHp coordinate
%   XYZ - vectors in XYZ coordinate
%
% NOTE:
%   This is used with [f_get_dw_orientation.m].
%
%
% Last modified
% 2010.08.06.
% 2010.10.22.
%   Take care of DWI case in which diffusion gradient orientation is one of
%   M, P, S.
%   Add slice_orientation to input.
% 2010.11.15.
%   Change order and sign of dw_ori_m based on object position and Matlab coordinate.
%   See note #4 page 17.
% 2010.12.09.
%   Add 'number_of_directions' for 'user defined' diffusion gradient
%   orientation.
%   Change 'dw_m' from which it uses GR`diff[0,1,2]:str_factor to which it
%   uses gradient orientations shown in PAR/REC file.
%
% 2011.01.19.
%   Use Philips coordinate relationship,
%       XYZ -> LPH <- L'P'H' <- NWV <- MPS
%   here, LPH is diff orientation in PAR file and NWV is actual image-space
%   coordinate of diff orientation. Then LPH (actually P,H,L order in PAR
%   file) must be transformed to NWV, and NWV must be re-ordered for CRP
%   (column, row, page) coordinate for tensor calculation in Matlab. When CRP
%   diff orientation is used, DW images must be in default display
%   orientation. This is shown in E4p40. This is tested in
%   [test_diffusion_gradient_coordinate.m].
% 2011.01.24.
%   This function [f_coord_tform.m] is generated from
%   [f_get_dw_orientation.m].



%% Check input
if size(PHL,1)~=3
    PHL = PHL';
end



%% Coordinate transform

% List of Tprep
Tprep_par = eye(3);
Tprep_per = [0 -1 0; 1 0 0; 0 0 1];


% List of Tfsd
Tfsd_m = [-1 0 0; 0 1 0; 0 0 1];
Tfsd_p = [1 0 0; 0 -1 0; 0 0 1];
Tfsd_s = [1 0 0; 0 1 0; 0 0 -1];


% List of Tsom
Tsom_tra = [0 -1 0; -1 0 0; 0 0 1];
Tsom_cor = [0 -1 0; 0 0 1; 1 0 0];
Tsom_sag = [0 0 -1; 0 -1 0; 1 0 0];


% Tprep, Tfsd, Tsom
if strcmpi(slice_orientation,'transverse')
    Tsom = Tsom_tra;
    
    switch fold_over_dir
        case 'AP'
            prep_s = 'per';
            switch fat_shift_dir
                case 'A'
                    fsd_s = 'm';
                case 'P'
                    fsd_s = 'p';
                otherwise
                    error('f_get_dw_orientation:main','Unknown fat_shift_dir')
            end
        case 'RL'
            prep_s = 'par';
            switch fat_shift_dir
                case 'R'
                    fsd_s = 'p';
                case 'L'
                    fsd_s = 'm';
                otherwise
                    error('f_get_dw_orientation:main','Unknown fat_shift_dir')
            end
        otherwise
            error('f_get_dw_orientation:main','Unknown fold_over_dir')
    end
    
elseif strcmpi(slice_orientation,'coronal')
    Tsom = Tsom_cor;
    
    switch fold_over_dir
        case 'RL'
            prep_s = 'par';
            switch fat_shift_dir
                case 'R'
                    fsd_s = 'p';
                case 'L'
                    fsd_s = 'm';
                otherwise
                    error('f_get_dw_orientation:main','Unknown fat_shift_dir')
            end
        case 'FH'
            prep_s = 'per';
            switch fat_shift_dir
                case 'F'
                    fsd_s = 'p';
                case 'H'
                    fsd_s = 'm';
                otherwise
                    error('f_get_dw_orientation:main','Unknown fat_shift_dir')
            end
        otherwise
            error('f_get_dw_orientation:main','Unknown fold_over_dir')
    end
    
elseif strcmpi(slice_orientation,'sagittal')
    Tsom = Tsom_sag;
    
    switch fold_over_dir
        case 'AP'
            prep_s = 'par';
            switch fat_shift_dir
                case 'A'
                    fsd_s = 'p';
                case 'P'
                    fsd_s = 'm';
                otherwise
                    error('f_get_dw_orientation:main','Unknown fat_shift_dir')
            end
        case 'FH'
            prep_s = 'per';
            switch fat_shift_dir
                case 'F'
                    fsd_s = 'p';
                case 'H'
                    fsd_s = 'm';
                otherwise
                    error('f_get_dw_orientation:main','Unknown fat_shift_dir')
            end
        otherwise
            error('f_get_dw_orientation:main','Unknown fold_over_dir')
    end
    
else
    error('f_get_dw_orientation:main','UNknown fold_over_dir')
end

if strcmpi(prep_s,'per')
    Tprep = Tprep_per;
else
    Tprep = Tprep_par;
end

if strcmpi(fsd_s,'m')
    Tfsd = Tfsd_m;
elseif strcmpi(fsd_s,'p')
    Tfsd = Tfsd_p;
else
    error('f_get_dw_orientation:main','Unknown fsd_s')
end


% Txyz
Txyz = Tprep*Tfsd;


% Trl, Tap, Tfh
Trl = @(rl)[1 0 0; ...
           0 cosd(rl) -sind(rl); ...
           0 sind(rl) cosd(rl)];
Tap = @(ap)[cosd(ap) 0 sind(ap); ...
           0 1 0; ...
           -sind(ap) 0 cosd(ap)];
Tfh = @(fh)[cosd(fh) -sind(fh) 0; ...
           sind(fh) cosd(fh) 0; ...
           0 0 1];


% Tang
AP = Ang.AP;   % degrees
RL = Ang.RL ;
FH = Ang.FH;

Tang = Trl(RL)*Tap(AP)*Tfh(FH);


% List of Tpo
Tpo_supine = eye(3);
Tpo_prone = [-1 0 0; 0 -1 0; 0 0 1];
Tpo_rightdecubitus = [0 -1 0; 1 0 0; 0 0 1];
Tpo_leftdecubitus = [0 1 0; -1 0 0; 0 0 1];


% List of Tpp
Tpp_headfirst = [0 1 0; -1 0 0; 0 0 -1];
Tpp_feetfirst = [0 -1 0; -1 0 0; 0 0 1];


% Tpom
Tpo_supine = eye(3);
Tpo_prone = [-1 0 0; 0 -1 0; 0 0 1];
Tpo_rightdecubitus = [0 -1 0; 1 0 0; 0 0 1];
Tpo_leftdecubitus = [0 1 0; -1 0 0; 0 0 1];
Tpp_feetfirst = [0 -1 0; -1 0 0; 0 0 1];
Tpp_headfirst = [0 1 0; -1 0 0; 0 0 -1];

switch patient_orientation
    case 'supine'
        Tpo = Tpo_supine;
    case 'prone'
        Tpo = Tpo_prone;
    case 'rightdecubitus'
        Tpo = Tpo_rightdecubitus;
    case 'leftdecubitus'
        Tpo = Tpo_leftdecubitus;
    otherwise
        error('f_get_dw_orientation:main','Unknown patient_orientation')
end

switch patient_position
    case 'feetfirst'
        Tpp = Tpp_feetfirst;
    case 'headfirst'
        Tpp = Tpp_headfirst;
    otherwise
        error('f_get_dw_orientation:main','Unknown patient_position')
end

Tpom = Tpo*Tpp;



%% Coordinate transformation
% There are five coordinate system in the scanner,
%   XYZ, LPH, L'P'H'(LPHp), NWV, MPS
% and in Matlab,
%   CRP (column,row,page)
LPH = PHL([3,1,2],:);
XYZ = inv(Tpom)*LPH;
LPHp = inv(Tang)*LPH;
NWV = inv(Tsom)*LPHp;
MPS = inv(Txyz)*NWV;



%% Take NWV as image-frame coordinate (default display orientation) and re-order NWV to CRP
% CRP is for default display orientation, see E4p44.

% Generate transform matrix from NWV to CRP.
Tcrp = [0 -1 0; -1 0 0; 0 0 1];

% Transform NWV to CRP.
CRP = Tcrp*NWV;



%% Final diffusion gradient orientation
CRP = CRP';
NWV = NWV';
MPS = MPS';
LPH = LPH';
LPHp = LPHp';
XYZ = XYZ';

% switch nargout
%     case 1
%         varargout{1} = CRP';
%     case 2
%         varargout{1} = CRP';
%         varargout{2} = NWV';
%     case 3
%         varargout{1} = CRP';
%         varargout{2} = NWV';
%         varargout{3} = MPS';
%     case 4
%         varargout{1} = CRP';
%         varargout{2} = NWV';
%         varargout{3} = MPS';
%         varargout{4} = LPH';
%     case 5
%         varargout{1} = CRP';
%         varargout{2} = NWV';
%         varargout{3} = MPS';
%         varargout{4} = LPH';
%         varargout{5} = LPHp';
%     case 6
%         varargout{1} = CRP';
%         varargout{2} = NWV';
%         varargout{3} = MPS';
%         varargout{4} = LPH';
%         varargout{5} = LPHp';
%         varargout{6} = XYZ';
%     otherwise
%         error('f_coord_tform:main','Unknown number of output')
% end



%% END





