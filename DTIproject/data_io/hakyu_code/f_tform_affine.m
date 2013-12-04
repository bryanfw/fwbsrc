
function T = f_tform_affine(data,p,method_s)
%[f_tform_affine] generates 12 parameter affine transformation matrix upto
%12 parameters.
%
% Usage:
%   T = f_tform_affine(data,p,method_s)
%
% Input:
%   data -
%       Input image
%   p -
%       Transformation parameters are, in general,
%       [trans(x,y,z),rot(x,y,z),scale(x,y,z),shear(x,y,z)]
%           trans -
%               positive number translates positive x,y axes in Matlab
%               image coordinate, in voxel
%           rot -
%               positive number rotates based on right-hand rule in Matlab
%               image coordinate, in degree
%           scale -
%               positive number magnifies image
%           shear -
%               positive number stretches top (x) and right (y) part of
%               image along positive x,y axes in Matlab image coordinate
%       but the number of p must be variable based on the dimension of data
%       and method_s.
%   method_s -
%       Method for co-registration,
%           'rigid'     rigid-body transformation
%           'affine'    affine transformation
%
% Output:
%   T -
%       Transform array
%
%
% Last modified
% 2011.01.19.
%   This function is generated from [tformAffine.m].
% 2011.01.20.
%   Referenced [spm_matrix.m].
% 2011.01.21.
%   Modified for 2-D and 3-D.
%
% Ha-Kyu



%% Check input
DIM = length(size(data));
if DIM > 3 || DIM < 2
    error('f_tform_affine:main','Image data must be 2D or 3D.')
end
[ny,nx,nz] = size(data); % 2-D or 3-D
if size(p,1) <= size(p,2)
    p = p';
end


% Data dimension, method and length of affine parameters.
if DIM==3
    if strcmpi(method_s,'rigid')
        if length(p) ~= 6
            error('f_tform_affine:main', ...
                '%d-D data, %s co-registration requires 6 parameters', ...
                DIM,upper(method_s))
        end
    elseif strcmpi(method_s,'affine')
        if length(p) ~= 12
            error('f_tform_affine:main', ...
                '%d-D data, %s co-registration requires 12 parameters', ...
                DIM,upper(method_s))
        end
    else
        error('f_tform_affine:coreg:main','Unknown method_s')
    end
elseif DIM==2
    if strcmpi(method_s,'rigid')
        if length(p) ~= 3
            error('f_tform_affine:main', ...
                '%d-D data, %s co-registration requires 3 parameters', ...
                DIM,upper(method_s))
        end
    elseif strcmpi(method_s,'affine')
        if length(p) ~= 7
            error('f_tform_affine:main', ...
                '%d-D data, %s co-registration requires 7 parameters', ...
                DIM,upper(method_s))
        end
    else
        error('f_tform_affine:coreg:main','Unknown method_s')
    end
else
    error('f_tform_affine:coreg:main','Dimension of data must be 2-D or 3-D')
end




%% Make generalized 12 parameter for all 2-D or 3-D transform
if DIM==3
    if strcmpi(method_s,'rigid') % 6 parameter
        p_temp = [p(1) p(2) p(3), p(4) p(5) p(6), 1 1 1, 0 0 0]';
        p = p_temp;
    elseif strcmpi(method_s,'affine') % 12 parameter
        % p = p;
        % Do nothing, it's already 12 parameter.
    end        
elseif DIM==2
    if strcmpi(method_s,'rigid') % 3 parameter
        p_temp = [p(1) p(2) 0, 0 0 p(3), 1 1 1, 0 0 0]';
        p = p_temp;
    elseif strcmpi(method_s,'affine') % 7 parameter
        p_temp = [p(1) p(2) 0, 0 0 p(3), p(4) p(5) 1, p(6) p(7) 0]';
        p = p_temp;
    end
end


% Check illegal scaling if scaling is less than or equal to 0.
% if DIM==3 && any(p(7:9)<=0)
%     error('f_tformAffine:main','Scaling parameter must be greater than zero.')
% end
% if DIM==2 && any(p(7:8)<=0)
%     error('f_tformAffine:main','Scaling parameter must be greater than zero.')
% end



%% Make transform array
%* It is in homogeneous coordinate system.
%* Positive input values for translation in p moves positive X and Y direction.
%* Rotations are all in right-hand coordinate system.
%* Scaling is proportional to the input values for scaling, p(7)=2, then 2 times in X direction.

if DIM==3
    % Translation.
    Mt = [1 0 0 -p(1); 0 1 0 -p(2); 0 0 1 -p(3); 0 0 0 1];
    
    % Rotation.
    Mrx = [1 0 0 0; 0 cosd(p(4)) sind(p(4)) 0; 0 -sind(p(4)) cosd(p(4)) 0; 0 0 0 1];    % X roatation
    Mry = [cosd(p(5)) 0 -sind(p(5)) 0; 0 1 0 0; sind(p(5)) 0 cosd(p(5)) 0; 0 0 0 1];    % Y roatation
    Mrz = [cosd(p(6)) sind(p(6)) 0 0; -sind(p(6)) cosd(p(6)) 0 0; 0 0 1 0; 0 0 0 1];    % Z roatation
    
    % Scaling.
    Msc = [1/p(7) 0 0 0; 0 1/p(8) 0 0; 0 0 1/p(9) 0; 0 0 0 1];
    
    % Shearing.
    %* 1.
    %Mshx = [1 p(10) 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    %Mshy = [1 0 0 0; p(11) 1 0 0; 0 0 1 0; 0 0 0 1];
    %Mshzx = [1 p(10) 0 0; 0 1 0 0; 0 p(12) 1 0; 0 0 0 1];
    %Mshzy = [1 0 0 0; p(11) 1 0 0; p(12) 0 1 0; 0 0 0 1];
    %Msh = Mshzy*Mshzx*Mshy*Mshx;
    
    %* 2.
    %Mshxy = [1 0 p(10) 0; 0 1 p(11) 0; 0 0 1 0; 0 0 0 1];
    %Msh = Mshxy*Mshy*Mshx;
    
    %* 3.
    %Msh = [1 0 p(10) 0; 0 1 p(11) 0; 0 0 1 0; 0 0 0 1];
    
    %* 4. Use this as [spm_matrix.m].
    Msh = [1 p(10) p(11) 0; 0 1 p(12) 0; 0 0 1 0; 0 0 0 1];
    
    
    % Translation - forward and backward.
    % Mtf = [1 0 0 -floor(nx/2); 0 1 0 -floor(ny/2); 0 0 1 -floor(nz/2); 0 0 0 1];  % forward translation
    % Mtb = [1 0 0 floor(nx/2); 0 1 0 floor(ny/2); 0 0 1 floor(nz/2); 0 0 0 1]; % backward translation
    Mtf = [1 0 0 -(nx/2); 0 1 0 -(ny/2); 0 0 1 -(nz/2); 0 0 0 1]; % forward translation
    Mtb = [1 0 0 (nx/2); 0 1 0 (ny/2); 0 0 1 (nz/2); 0 0 0 1]; % backward translation
    
else
    % Translation.
    Mt = [1 0 -p(1); 0 1 -p(2); 0 0 1];
    
    % Rotation.
    Mrx = eye(3);
    Mry = eye(3);
    Mrz = [cosd(p(6)) sind(p(6)) 0; -sind(p(6)) cosd(p(6)) 0; 0 0 1];    % Z roatation
    
    % Scaling.
    Msc = [1/p(7) 0 0; 0 1/p(8) 0; 0 0 1];
    
    % Shearing.
    Msh = [1 p(10) 0; p(11) 1 0; 0 0 1];
    
    % Translation - forward and backward.
    Mtf = [1 0 -(nx/2); 0 1 -(ny/2); 0 0 1]; % forward translation
    Mtb = [1 0  (nx/2); 0 1  (ny/2); 0 0 1]; % backward translation
    
end


% Final transform matrix.
T1 = Mt*Mrx*Mry*Mrz*Msc*Msh;
T = Mtb*T1*Mtf;
clear  T1




%% END




