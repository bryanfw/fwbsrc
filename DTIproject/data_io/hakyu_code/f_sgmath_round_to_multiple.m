
function val_out = f_sgmath_round_to_multiple(val_in, step_size)
%[f_sgmath_round_to_multiple] find round up value to the multiple of
%step_size.
%
% USAGE:
%   val_out = f_sgmath_round_to_multiple(val_in, step_size)
%
% NOTE:
%   This function originates from [SGMATH_round_to_multiple.m].
%
%
% Last modified
% 2010.09.01.
%   Modified based on 3T R2.6.3.4 SWID 118,
%   ugeo0_aspect_ratio-[mpugeo__g.c].



%% Calculate the nearest multiple of step_size from input
if mod(val_in,step_size)==0
    val_out = val_in;
else
    val_mod = mod(val_in,step_size);
    if (step_size-val_mod) > (step_size/2)
        val_out = val_in - val_mod;
    else
        val_out = val_in - val_mod + step_size;
    end
end
if mod(val_out,step_size)~=0
    error('Wrong procedure.')
end

return








    
    