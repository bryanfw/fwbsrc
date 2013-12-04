
%[s_get_GENparams] get general parameters, GENparams.
%
% USAGE:
%   s_get_GENparams
%
% OUTPUT:
%   GENparams: A struct with fields as,
%       B0, in Gauss
%       gamma, in rad/(mT*s)
%
%
% Last modified
% 2010.08.06.
% 2010.09.02.
%   Take B0.
% 2010.09.06.
%   Take coilID from EX_GEO_cur_stack_coil_id defined in batch file.
% 2010.09.16.
%   Don't clear variables. This will be cleared later.
% 2011.01.19.
%   Add 'patient_orientation' and 'patient_position'.
%
% Ha-Kyu



%% Define GENparams

% Clear.
clear  GENparams

% GENparams.
GENparams.B0 = main_field_strength; % Gauss
GENparams.gamma = 42577.46778*2*pi; % rad/(mT*s)
GENparams.coilID = coilType; % EX_GEO_cur_stack_coil_id
GENparams.os = computer; % PC or not (Linux, Mac or Solaris etc.)
GENparams.patient_orientation = patient_orientation;
GENparams.patient_position = patient_position;

% Clear.
% clear  main_field_strength








