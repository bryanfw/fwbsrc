
function varargout = f_get_ENCima_v2(scan_resol, recon_resol, sense_factor, isEPI, ...
    epi_factors, slice_thickness, fov, UGN1_ACQ_reduced_acq_time, aspect_ratio,DWIparams)
% [f_get_ENCima_v2] calculates ENC`ima (almost fully), LCAima, STACKima, ACQREC
% objects.
%
% USAGE
%   f_get_ENCima_v2(scan_resol, recon_resol, sense_factor, isEPI, epi_factors, 
%                 slice_thickness, fov, UGN1_ACQ_reduced_acq_time, aspect_ratio,DWIparams)
%
% INPUT
%   scan_resol                  :   scan resolution
%   recon_resol                 :   recon resolution
%   sense_factor                :   sense factor in P_ORI (user input value)
%   isEPI                       :   true or false for EPI sequence
%   epi_factors                 :   EPI factor
%   slice_thickness             :   slice thickness, mm
%   fov                         :   Field of view, mm
%   UGN1_ACQ_reduced_acq_time   :   `UGN1_ACQ_reduced_acq_time = 
%                                   `EX_ACQ_reduced_acq_time, [mpuacq__g.c]
%   aspect_ratio                :   STACK`ima:aspect_ratio. This is not 1 when
%                                   rectangular FOV is used.
%
%
% Last modified
%   2009.02.11
%   2009.02.16
%   2009.02.24.     Add varargout.
%   2009.05.12.     Add output for ACQ_matrix_MxP.
%   2009.09.02.     Modify actual SENSE FACTOR.
%                   Add slice thickness and fov.
%   2009.10.22.     Make all VARARGIN as input arguments.
%                   Take `UGN1_ACQ_reduced_acq_time as an input because this changes
%                   with voxel size, especially at 3T.
%   2010.----.      Add LCA`ima:[echo]:gamma_fovs.
%   2010.07.09.     Take STACK'ima:aspect_ratio as an input for rectangular FOV.
%   2010.07.12.     Change ACQREC.number_of_shots from uint8 to double.
% 2010.08.06.
%   Function name changed to from [get_ENCima_v3.m] to
%   [f_get_ENCima.m] with prefix 'f' to show that this is function.
% 2010.09.01.
%   SGMATH_round_to_multiple() is modified and generated as external
%   function, [f_sgmath_round_to_multiple.m].
% 2012.07.23.
%   This version 2 is generated from [f_get_ENCima.m].
%
%
% See also,     get_ENCima, get_ENCima_v2, f_get_ENCima
%
% Ha-Kyu



%% Parse input    
    
    
%% Variable objects
UGN1_ACQ_scan_resol = scan_resol;
UGN1_PROC_recon_resol = recon_resol;
UGN1_GEO_sense_m_red_factor = 1;
UGN1_GEO_sense_p_red_factor = sense_factor;      % use input sense factor
UGN1_ACQ_measurements = 1.0;
UGN1_ACQ_act_measurements = 1.0;
% UGN1_ACQ_reduced_acq_time = 100;      % at 7T it always is 100
% UGN1_ACQ_reduced_acq_time = 81.33;    % at 3T it changes with voxel size
ICSC_sense_basic_ovs_margin = 0.5;
ICSC_sense_gibbs_ovs_margin = 2.0;


%% Constant objects
M_ORI = 1;
P_ORI = 2;
S_ORI = 3;
step_size = 4;


%% Other objects
init_const = nan;
init_v = ones(1,3)*nan;
STACKima = struct('aspect_ratio',init_const);
LCAima = struct('pre_sense_samples',init_v,'samples',init_v);
ENCima = struct('oversample_resolutions',init_v, ...
    'oversample_factors',init_v,'scan_resolutions',init_v, ...
    'profiles',init_v,'interp_factors',init_v,'ft_lengths',init_v, ...
    'recon_resolutions',init_v,'epi_factors',init_const,...
    'partial_matrix_factors',init_v,'sense_factors',init_v, ...
    'max_encoding_numbers',init_v,'min_encoding_numbers',init_v);

% Initialize
%* This is the same as [RFOV (%)] in [GEOMETRY] tab in ExamCard (not in console).
% STACKima.aspect_ratio = 1.0;
% STACKima.aspect_ratio = 0.8;
STACKima.aspect_ratio = aspect_ratio;


%% Read .LIST file for some variables
% After 2ch multix upgrade to 3Tb, there is an occasion for ovs_fac(M_ORI)
% is 4, not usual value of 2 in [scan20120713_r3099__Smith_207380_3T.m].
scan_params = f_get_params_from_listfile( ...
    fullfile(pwd,[DWIparams.filename,'.list']));
pack


%% Calculate ENC`ima, LCA`ima and other scan parameters (say, ACQREC)

%------------------------
% sense_factors: this will be modified later.
ENCima.sense_factors(M_ORI) = scan_params.X_direction_SENSE_factor;%1.0;
ENCima.sense_factors(P_ORI) = sense_factor;
ENCima.sense_factors(S_ORI) = 1.0;

%------------------------
% epi_factors
ENCima.epi_factors = epi_factors;

%------------------------
% recon_resolutions
recon_res = init_v;
recon_res(M_ORI) = UGN1_PROC_recon_resol;
recon_res(P_ORI) = max(1, f_sgmath_round_to_multiple( ...
    UGN1_PROC_recon_resol * STACKima.aspect_ratio, step_size) );
recon_res(S_ORI) = 1;

%%%
ENCima.recon_resolutions = recon_res;
%%%

%------------------------
% oversample_factors
sense_ovs_factor = init_v;
sense_ovs_factor(M_ORI) = scan_params.kx_oversample_factor;%2.0;
hs_margin = 0.0;
margin_pixels = uint16( round(ICSC_sense_basic_ovs_margin + ...
    max(hs_margin, ICSC_sense_gibbs_ovs_margin)) );
nominal_pixels = uint16( round(UGN1_ACQ_scan_resol * ...
    STACKima.aspect_ratio) );
sense_ovs_factor(P_ORI) = min( double((nominal_pixels+2*margin_pixels)) / ...
    double(nominal_pixels), UGN1_GEO_sense_p_red_factor );
sense_ovs_factor(S_ORI) = 1;

%%%
LCAima.oversample_factors(M_ORI) = sense_ovs_factor(M_ORI);
LCAima.oversample_factors(P_ORI) = max(sense_ovs_factor(P_ORI), ...
    UGN1_ACQ_measurements/UGN1_ACQ_act_measurements);
LCAima.oversample_factors(S_ORI) = 1.0;
ENCima.oversample_factors = LCAima.oversample_factors;
%%%

%------------------------
% partial_matrix_factors
partial_matrix_factors = ones(1,3);

%%%
ENCima.partial_matrix_factors(M_ORI) = partial_matrix_factors(M_ORI);
ENCima.partial_matrix_factors(P_ORI) = partial_matrix_factors(P_ORI);
ENCima.partial_matrix_factors(S_ORI) = partial_matrix_factors(S_ORI);
%%%

%------------------------
% samples
pre_sense_samples = init_v;
samples = init_v;

pre_sense_samples(M_ORI) = UGN1_ACQ_scan_resol * ...
    ENCima.oversample_factors(M_ORI);
pre_sense_samples(P_ORI) = UGN1_ACQ_scan_resol * STACKima.aspect_ratio * ...
    ENCima.oversample_factors(P_ORI);
pre_sense_samples(S_ORI) = 1;

samples(M_ORI) = min(pre_sense_samples(M_ORI), ...
    ceil(pre_sense_samples(M_ORI)/UGN1_GEO_sense_m_red_factor));
samples(P_ORI) = min(pre_sense_samples(P_ORI), ...
    ceil(pre_sense_samples(P_ORI)/UGN1_GEO_sense_p_red_factor));
samples(S_ORI) = pre_sense_samples(S_ORI)*LCAima.oversample_factors(S_ORI);

%%%
LCAima.pre_sense_samples = pre_sense_samples;
LCAima.samples = samples;
% LCAima.pre_sense_samples = LCAima.samples;
%%%

%------------------------
% sense factors
LCAima.sense_factors = LCAima.pre_sense_samples./LCAima.samples;
ENCima.sense_factors = LCAima.sense_factors;

%------------------------
% oversample_resolutions
ENCima.oversample_resolutions(M_ORI) = LCAima.samples(M_ORI);
ENCima.oversample_resolutions(P_ORI) = round(LCAima.samples(P_ORI) * ...
    UGN1_ACQ_reduced_acq_time/100);
ENCima.oversample_resolutions(S_ORI) = LCAima.samples(S_ORI);
ENCima.profiles = ENCima.oversample_resolutions;

if isEPI
    py_factor = uint16(ENCima.epi_factors);
    py_shots = floor(ENCima.profiles(P_ORI)/double(py_factor));
    ovs_res = min(round(py_factor*py_shots / ...
        ENCima.partial_matrix_factors(P_ORI)), ...
        ENCima.oversample_resolutions(P_ORI));
    ENCima.oversample_resolutions(P_ORI) = ovs_res;
    ENCima.profiles(P_ORI) = ENCima.oversample_resolutions(P_ORI);
end

%------------------------
% max(min)_encoding_numbers
for ind_ORI = [M_ORI, P_ORI, S_ORI]    
    ENCima.max_encoding_numbers(ind_ORI) = ...
        floor( (ENCima.oversample_resolutions(ind_ORI)-1)/2 );
    ENCima.min_encoding_numbers(ind_ORI) = ...
        -( ENCima.oversample_resolutions(ind_ORI) - ...
        ENCima.max_encoding_numbers(ind_ORI) - 1 );
end
        
%------------------------
% scan_resolutions
for ind_ORI = [M_ORI, P_ORI, S_ORI]    
    ENCima.scan_resolutions(ind_ORI) = ...
        ENCima.oversample_resolutions(ind_ORI) / ...
        ENCima.oversample_factors(ind_ORI);
end

%------------------------
% interp_factors
for ind_ORI = [M_ORI, P_ORI, S_ORI]    
    ENCima.interp_factors(ind_ORI) = ...
        ENCima.recon_resolutions(ind_ORI) * ...
        ENCima.oversample_factors(ind_ORI) / ...
        ( ENCima.oversample_resolutions(ind_ORI) * ...
        ENCima.sense_factors(ind_ORI) );
end

%------------------------
% ft_lengths
for ind_ORI = [M_ORI, P_ORI, S_ORI]
    ENCima.ft_lengths(ind_ORI) = ENCima.recon_resolutions(ind_ORI) * ...
        ENCima.oversample_factors(ind_ORI) / ...
        ENCima.sense_factors(ind_ORI);
end

%------------------------
% ACQ matrix MxP
ACQ_matrix_MxP = floor( ENCima.scan_resolutions([M_ORI,P_ORI]) .* ...
    ENCima.sense_factors([M_ORI,P_ORI]) );

%------------------------
% ACQ voxel MPS
ACQ_voxel_MPS = fov./(ENCima.scan_resolutions(1:2).*ENCima.sense_factors(1:2));
ACQ_voxel_MPS(S_ORI) = slice_thickness;

%------------------------
% REC voxel MPS
%$REC_voxel_MPS = repmat(fov/recon_resol,1,2);
%$REC_voxel_MPS = fov./[recon_resol,recon_resol];
REC_voxel_MPS = ACQ_voxel_MPS./ENCima.interp_factors;
REC_voxel_MPS(S_ORI) = slice_thickness;

%------------------------
% Number of shots
% number_of_shots = uint8(ENCima.oversample_resolutions(P_ORI)/epi_factors);
number_of_shots = double(uint8(ENCima.oversample_resolutions(P_ORI)/epi_factors));

%------------------------
% ACQREC
ACQREC.ACQ_matrix_MxP = ACQ_matrix_MxP;
ACQREC.ACQ_voxel_MPS = ACQ_voxel_MPS;
ACQREC.REC_voxel_MPS = REC_voxel_MPS;
ACQREC.number_of_shots = number_of_shots;


%% Report
if nargout==0
    fprintf('\n\n')
    disp('----------------------------------------------------------------')    

    fprintf('ENC`ima\n')
    disp(ENCima)
    
    fprintf('\n')
    fprintf('LCA`ima\n')
    disp(LCAima)
    
    fprintf('\n')
    fprintf('STACK`ima\n')
    disp(STACKima)
    
    fprintf('\n')
    fprintf('ACQ matrix MxP\n')
    disp(ACQ_matrix_MxP)
    
    fprintf('\n')
    fprintf('ACQ voxel MPS\n')
    disp(round(100* ACQ_voxel_MPS )/100)     % round at 3rd digit under zero
    fprintf('    (%f, %f, %f)\n',ACQ_voxel_MPS)
    
    fprintf('\n')
    fprintf('REC voxel MPS\n')
    disp(REC_voxel_MPS)
    
    fprintf('\n')
    fprintf('Number of Shots\n')
    disp(number_of_shots)
    
    disp('----------------------------------------------------------------')
    fprintf('\n')
    
elseif nargout==1
    varargout{1} = ENCima;
elseif nargout==2
    varargout{1} = ENCima;
    varargout{2} = LCAima;
elseif nargout==3
    varargout{1} = ENCima;
    varargout{2} = LCAima;
    varargout{3} = STACKima;
elseif nargout==4
    varargout{1} = ENCima;
    varargout{2} = LCAima;
    varargout{3} = STACKima;
    varargout{4} = ACQREC;
else
    error('Too many nargout.')
end


% %% Subfunctions
% function val_out = SGMATH_round_to_multiple(val_in, step_size)
% % Find round up value to the multiple of step_size.
% 
% if mod(val_in,step_size)==0
%     val_out = val_in;
% else
%     %val_mod = mod(val_in,step_size);
%     %val_out = val_in - val_mod + step_size;
%     
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
    



%% END












