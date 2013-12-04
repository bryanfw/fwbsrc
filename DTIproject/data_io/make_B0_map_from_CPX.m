%% MAKE_B0_MAP_FROM_CPX     Create B0 Map in Hz from a vuSuperExport ZIP file that contains a Philips CPX file
%
% FILENAME_B0_MAP_NEW = MAKE_B0_MAP_FROM_CPX(FILENAME_VUSUPEREXPORT)
%
%   FILENAME_VUSUPEREXPORT is a string containing a file prefix or name of 
%   the vuSuperExport ZIP file, CPX file, PAR file or REC file, e.g. 
%   MYSERIES or MYSERIES.ZIP or MYSERIES.CPX or MYSERIES.PAR or
%   MYSERIES.REC
%
%   The ZIP file should contain a Philips CPX file and PARREC file pair for
%   the magnitude-only (incorrect) B0 Map.
%
%   Also, the delta TE time is not available in the PAR file and is
%   therefore taken from EX_ACQ_B0_map_delta_TE in the .PDF.XML file that 
%   should be present in the vuSuperExport ZIP file.
%
%   FILENAME_B0_MAP_NEW is a string containing the name of the newly
%   created PAR file, which will be named with the same prefix as the 
%   FILENAME_ZIP with "_B0_MAP_FROM_CPX" appended to the prefix, e.g.
%   MYSERIES_B0_MAP_FROM_CPX.PAR
%
%   If the function fails, FILENAME_B0_MAP_NEW will be an empty string.
%
% FILENAME_B0_MAP_NEW = MAKE_B0_MAP_FROM_CPX(FILENAME_VUSUPEREXPORT,MASK_FLAG)
%
%   If MASK_FLAG is non-zero, the B0 map will be masked by the mean
%   magnitude image after thresholding by Otsu's method. By default, the
%   B0 map is not masked.
%
%  See also: LOADPARREC, LOADCPX, WRITEPARREC, DICTPARREC, LOADPDFXML
%

%% Revision History
% * 2012.07.28    initial version - welcheb
function filename_B0_map_new = MAKE_B0_MAP_FROM_CPX(filename_VUSUPEREXPORT,mask_flag)

%% Initialize returned filename to an empty string
filename_B0_map_new = '';

%% Check number of input arguments
if nargin<1,
    error(sprintf('MAKE_B0_MAP_FROM_CPX requires at least 1 input argument.\n\nUsage: filename_B0_map_new = make_B0_map_from_CPX(filename_VUSUPEREXPORT);\n'));
    return;
end
if nargin<2,
    mask_flag = 0;
end

%% Parse the filename_VUSUPEREXPORT filename input
% It may be the ZIP filename, CPX filename, PAR filename, REC filename or just the filename prefix
% Instead of REGEXP, use REGEXPI which igores case
filename_VUSUPEREXPORT = regexprep(filename_VUSUPEREXPORT,'\\*','/');
filename_VUSUPEREXPORT = regexprep(filename_VUSUPEREXPORT,'/*','/');
toks = regexpi(filename_VUSUPEREXPORT,'^(.*?)(\.ZIP|\.CPX|\.PAR|\.REC)?$','tokens');
fileprefix_full = toks{1}{1};
[fileparentfolder, fileprefix_short_a, fileprefix_short_b] = fileparts(fileprefix_full);
if length(fileparentfolder)<1,
    fileparentfolder ='.';
end
fileprefix = sprintf('%s%s', fileprefix_short_a, fileprefix_short_b);
zipname_VUSUPEREXPORT = sprintf('%s.ZIP', fileprefix);
cpxname_VUSUPEREXPORT = sprintf('%s.CPX', fileprefix);
parname_VUSUPEREXPORT = sprintf('%s.PAR', fileprefix);
recname_VUSUPEREXPORT = sprintf('%s.REC', fileprefix);
pdfxmlname_VUSUPEREXPORT = sprintf('%s.PDF.XML', fileprefix);

zipname_VUSUPEREXPORT_fullpath = sprintf('%s/%s.ZIP', fileparentfolder, fileprefix);
cpxname_VUSUPEREXPORT_fullpath = sprintf('%s/%s.CPX', fileparentfolder, fileprefix);
parname_VUSUPEREXPORT_fullpath = sprintf('%s/%s.PAR', fileparentfolder, fileprefix);
recname_VUSUPEREXPORT_fullpath = sprintf('%s/%s.REC', fileparentfolder, fileprefix);
pdfxmlname_VUSUPEREXPORT_fullpath = sprintf('%s/%s.PDF.XML', fileparentfolder, fileprefix);

%% Now we can create parname_B0_map_new
parname_B0_map_new = sprintf('%s/%s_B0_MAP_FROM_CPX.PAR', fileparentfolder, fileprefix);

%% Make sure unzipped files are available

% CPX
if exist(cpxname_VUSUPEREXPORT_fullpath)~=2,
    system_command_str = sprintf('unzip "%s" "%s" -d "%s"', zipname_VUSUPEREXPORT_fullpath, cpxname_VUSUPEREXPORT, fileparentfolder);
    [status,result] = system(system_command_str);
    if exist(cpxname_VUSUPEREXPORT_fullpath)~=2,
        error(sprintf('Cannot unzip %s from %s', cpxname_VUSUPEREXPORT, zipname_VUSUPEREXPORT_fullpath));
    end
end

% PAR
if exist(parname_VUSUPEREXPORT_fullpath)~=2,
    system_command_str = sprintf('unzip "%s" "%s" -d "%s"', zipname_VUSUPEREXPORT_fullpath, parname_VUSUPEREXPORT, fileparentfolder);
    [status,result] = system(system_command_str);
    if exist(parname_VUSUPEREXPORT_fullpath)~=2,
        error(sprintf('Cannot unzip %s from %s', parname_VUSUPEREXPORT, zipname_VUSUPEREXPORT_fullpath));
    end
end

% REC
if exist(recname_VUSUPEREXPORT_fullpath)~=2,
    system_command_str = sprintf('unzip "%s" "%s" -d "%s"', zipname_VUSUPEREXPORT_fullpath, recname_VUSUPEREXPORT, fileparentfolder);
    [status,result] = system(system_command_str);
    if exist(recname_VUSUPEREXPORT_fullpath)~=2,
        error(sprintf('Cannot unzip %s from %s', recname_VUSUPEREXPORT, zipname_VUSUPEREXPORT_fullpath));
    end
end

% PDFXML
if exist(pdfxmlname_VUSUPEREXPORT_fullpath)~=2,
    system_command_str = sprintf('unzip "%s" "%s" -d "%s"', zipname_VUSUPEREXPORT_fullpath, pdfxmlname_VUSUPEREXPORT, fileparentfolder);
    [status,result] = system(system_command_str);
    if exist(pdfxmlname_VUSUPEREXPORT_fullpath)~=2,
        error(sprintf('Cannot unzip %s from %s', pdfxmlname_VUSUPEREXPORT, zipname_VUSUPEREXPORT_fullpath));
    end
end

%% Load PDF.XML
info_pdfxml = loadPDFXML(pdfxmlname_VUSUPEREXPORT_fullpath);

%% Load CPX
loadopts_cpx.ec = [0];
loadopts_cpx.dyn = [0];
loadopts_cpx.ph = [0];
loadopts_cpx.row = [0 1]; % different TEs stored as rows
loadopts_cpx.mix = [0];
loadopts_cpx.savememory = true;
loadopts_cpx.verbose = false;
[data_cpx,info_cpx] = loadCPX(cpxname_VUSUPEREXPORT_fullpath,loadopts_cpx);
data_cpx = squeeze(data_cpx);

%% Load Old PARREC
loadopts_parrec.ty = [0];
loadopts_parrec.seq = [2];
loadopts_parrec.scale = 'FP';
loadopts_parrec.savememory = true;
loadopts_parrec.verbose = false;
loadopts_parrec.reducesingletons = false;
[data_parrec,info_parrec] = loadPARREC(parname_VUSUPEREXPORT_fullpath,loadopts_parrec);

%% Calculate CPX data phase difference
if length(info_cpx.dims.coil)==1,
    % single channel data
    data_phase_difference = angle( data_cpx(:,:,:,2) .* conj(data_cpx(:,:,:,1)) );
    data_mean_magnitude = ( abs(data_cpx(:,:,:,1)) + abs(data_cpx(:,:,:,2)) )/2;
else
    % multi channel data
    nCoils = length(info_cpx.dims.coil);
    idx_coil = info_cpx.dims.coil+1; % coil indices in CPX start at 0
    
    % try to get noise estimates from info_pdfxml
    noise_sigma_real = info_pdfxml.PDF_PREP_PARS.GSR_RC_noise_real_arr;
    noise_sigma_imag = info_pdfxml.PDF_PREP_PARS.GSR_RC_noise_imaginary_arr;
    noise_variance = (noise_sigma_real.^2 + noise_sigma_imag.^2).^(0.5);
    
    % if not all noise estimates are non-zero, just set them all to 1
    if ~all(noise_variance(idx_coil)),
        noise_variance(idx_coil)=1;
    end

    % use method from Bernstein MRM 32:330-334 (1994), Equations and [13] and [14]
    tmpsum = complex(zeros(size(data_cpx(:,:,:,1,1))));
    for k=1:nCoils,
        tmpsum = tmpsum + squeeze( data_cpx(:,:,:,2,k) .* conj(data_cpx(:,:,:,1,k)) )./noise_variance(idx_coil(k));
    end
    data_phase_difference = angle(tmpsum);
    data_mean_magnitude = abs(tmpsum).^(0.5);
end

%% Flip each slice (1st dimension)
data_phase_difference = flipdim(data_phase_difference, 1);

%% Mask phase using Otsu's method on the mean magnitude image
if mask_flag,
    data_mean_magnitude = data_mean_magnitude - min(data_mean_magnitude(:));
    data_mean_magnitude = data_mean_magnitude / max(data_mean_magnitude(:));
    otsu_normalized_threshold = graythresh(data_mean_magnitude(:));
    data_mask_zero = data_mean_magnitude < otsu_normalized_threshold;
    data_mask_zero_idx = find(data_mask_zero(:));
    data_phase_difference(data_mask_zero_idx) = 0;
end

%% Convert radians to Hz
if isfield(info_pdfxml.PDF_EXAM_PARS,'EX_ACQ_B0_map_delta_TE'),
    delta_TE_sec = info_pdfxml.PDF_EXAM_PARS.EX_ACQ_B0_map_delta_TE / 1000; 
elseif isfield(info_pdfxml.PDF_EXAM_PARS,'EX_ACQ_delta_TE'),
    delta_TE_sec = info_pdfxml.PDF_EXAM_PARS.EX_ACQ_delta_TE / 1000;
else
    error('Cannot determine B0 Map delta_TE_sec from .PDF.XML file');
    return;
end
max_offresonance_Hz = (1/delta_TE_sec)/2;
data_Hz = max_offresonance_Hz * data_phase_difference/pi;

%% Modify PARREC data array
% Fiedmap is added in the last dimension
data_parrec_new = zeros([info_parrec.datasize(1:3) 2],'single');
data_parrec_new(:,:,:,1) = data_parrec;
data_parrec_new(:,:,:,2) = data_Hz;

%% Modify PARREC info structure
% B0 map data is ty 3 (phase), seq 4
% Scale to have minimum value of -max_offresonance_Hz
% writePARREC uses info.table_row_index_array and info.imgdef components for writing
info_parrec_new = info_parrec;
nSlices = info_parrec_new.datasize(3);
nTableRows_Old = size(info_parrec.table,1);
idx_new_tablerows = [1:nSlices]+nTableRows_Old;
% Append new row numbers to the table_row_index_array
info_parrec_new.table_row_index_array = [info_parrec_new.table_row_index_array(:) ; idx_new_tablerows(:) ];

% Loop through all imgdef fields and append to the vals arrays
idx_copy_rows = info_parrec_new.table_row_index_array(1:nSlices);
imgdef_fieldnames = fieldnames(info_parrec_new.imgdef);
for k=1:length(imgdef_fieldnames),
    old_vals = info_parrec_new.imgdef.(imgdef_fieldnames{k}).vals;
    copy_vals = old_vals(idx_copy_rows(:),:);
    info_parrec_new.imgdef.(imgdef_fieldnames{k}).vals = [ old_vals ; copy_vals];
end

% Number of images should double
info_parrec_new.n_data_imgs = 2 * info_parrec_new.n_data_imgs;

% Make specific changes to the table rows refering to the field map
info_parrec_new.imgdef.image_type_mr.vals(idx_new_tablerows) = 3;
info_parrec_new.imgdef.scanning_sequence.vals(idx_new_tablerows) = 4;
info_parrec_new.imgdef.rescale_intercept.vals(idx_new_tablerows) = -max_offresonance_Hz;
info_parrec_new.imgdef.rescale_slope.vals(idx_new_tablerows) = (2*max_offresonance_Hz)/4095;
info_parrec_new.imgdef.scale_slope.vals(idx_new_tablerows) = 4095/(2*max_offresonance_Hz);
info_parrec_new.imgdef.window_center.vals(idx_new_tablerows) = 0;
info_parrec_new.imgdef.window_width.vals(idx_new_tablerows) = 2 * max_offresonance_Hz;

%% Write PARREC
if writePARREC(parname_B0_map_new, data_parrec_new, info_parrec_new, '4.2'),
    filename_B0_map_new = parname_B0_map_new;
else
    error(sprintf('writePARREC not successful'));
end