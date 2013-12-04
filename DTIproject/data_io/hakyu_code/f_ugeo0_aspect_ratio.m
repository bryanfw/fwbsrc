
function aspect_ratio = f_ugeo0_aspect_ratio(EX_ACQ_scan_resol,EX_GEO_rect_fov_perc)
%[f_ugeo0_aspect_ratio] calculates aspect ratio of FOV. This aspect ratio
%will be used to STACK`ima.aspect_ratio which can be acquired in examcard
%of 7T not in 3T.
%
% USAGE
%   aspect_ratio = f_ugeo0_aspect_ratio(EX_ACQ_scan_resol,EX_GEO_rect_fov_perc)
%
%
% Modified
% 2010.09.01.
%   This function is taken from ugeo0_aspect_ratio-[mpugeo__g.c] of 3T
%   R2.6.3.4_SWID_118.
%   SGMATH_round_to_multiple() is generated as an external function,
%   [f_sgmath_round_to_multiple.m].
% HKJ



%% Calculate aspect ratio
fov_perc = EX_GEO_rect_fov_perc/100;
step_size = 4;
%y_resol = SGMATH_round_to_multiple(EX_ACQ_scan_resol*fov_perc,step_size);
y_resol = f_sgmath_round_to_multiple(EX_ACQ_scan_resol*fov_perc,step_size);
y_resol = max(step_size,y_resol);
aspect_ratio = y_resol / double(EX_ACQ_scan_resol);



% %% Subfunctions
% function val_out = SGMATH_round_to_multiple(val_in, step_size)
% % Find round up value to the multiple of step_size.
% 
% if mod(val_in,step_size)==0
%     val_out = val_in;
% else
%     val_mod = mod(val_in,step_size);
%     if (step_size-val_mod) > (step_size/2)
%         val_out = val_in - val_mod;
%     else
%         val_out = val_in - val_mod + step_size;
%     end
% end
% if mod(val_out,step_size)~=0
%     error('Wrong procedure.')
% end
% 
% return
    





