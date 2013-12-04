
function ky = f_get_nav_ky(Eimg,Enav,Ns)
%[f_get_nav_ky] get navigator-echo ky position in image-echo space.
%
% Usage:
%   ky = f_get_nav_ky(Eimg,Enav,Ns)
%
% Input:
%   Eimg - image-echo epi-factor
%   Enav - navigator-echo epi-factor
%   Ns - number of shots
%
% Output:
%   ky - ky index of navigator-echo in image-echo space (Ns*Eimg lines)
%
%
% Last modified
% 2011.06.21.



%% Main
ky = floor(Eimg*Ns/2+1-floor(Enav/2)):floor(Eimg*Ns/2+1+ceil(Enav/2)-1);

