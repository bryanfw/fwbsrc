%
% V3
%
% Provide a version specification for reading the general information
% lines provided in a Philips NV .PAR file written by the research image
% export tool V3.
%
%
% VERSION SPECIFICATION
%
% Philips .PAR files are versioned. The version number of the 
% target file is contained in a comment line such as:
%
% # CLINICAL TRYOUT      Research image export tool  V3 
%
% Indicating which version of the image export tool that wrote
% the target file. Different versions differ in the number of
% general information lines they contain and the length and 
% content of the tag lines that describe the binary image data 
% in the corresponding .REC file. To make this read_parfile function 
% more flexible, format differences between different versions are
% encoded as a MATLAB cell array in separate M-FILES. The
% experienced user may create his own M-FILE specifications. 
%
% The cell array has the format:
%
% { line_key_string   value_type_string    field_name_string; ...
%   line_key_string   value_type_string    field_name_string; ...
%   line_key_string   value_type_string    field_name_string; ...
%   line_key_string   value_type_string    field_name_string; ...
% };
%
% Where the line_key_string matches a general information line tag such 
% as:
%
% '.    Examination name                   :'     
%
% and the value_type_string describes the type of data on the general
% information line that follows the colon delimeter. Valid 
% value_type_strings are:
%
% 'char  vector'         format string data as an cell array of words.   
% 'char  scaler'         format string data as a single string.    
% 'float vector'         format numeric data as an array of numeric values.
% 'float scaler'         format numeric data as a single value.   
%
% The value_type_strings determine how the information on the rest of
% the specified general information line is formated in the parfile
% structure.
%
% The field_name_strings can be any suitable memnonic for the extracted
% data. The parfile will contain fields with names specified by these
% strings.
%
% $Id: V3.m,v 1.1 2007/11/26 21:52:17 dfoxall Exp $
%

geninfo_spec = { ...                                  
'.    Patient name                       :'    'char  scalar'    'patient';    ...   
'.    Examination name                   :'    'char  vector'    'exam_name';   ... 
'.    Protocol name                      :'    'char  vector'    'protocol';   ... 
'.    Examination date/time              :'    'char  vector'    'exam_date';  ...
'.    Acquisition nr                     :'    'int   scalar'    'acq_nr';    ...
'.    Reconstruction nr                  :'    'int   scalar'    'recon_nr';  ...
'.    Scan Duration [sec]                :'    'float scalar'    'scan_dur';        ...
'.    Max. number of cardiac phases      :'    'int   scalar'    'max_card_phases'; ...
'.    Max. number of echoes              :'    'int   scalar'    'max_echoes'; ...
'.    Max. number of slices/locations    :'    'int   scalar'    'max_slices'; ... 
'.    Max. number of dynamics            :'    'int   scalar'    'max_dynamics'; ... 
'.    Max. number of mixes               :'    'int   scalar'    'max_mixes'; ... 
'.    Image pixel size [8 or 16 bits]    :'    'int   scalar'    'pixel_bits'; ... 
'.    Technique                          :'    'char  scalar'    'technique'; ...  
'.    Scan mode                          :'    'char  scalar'    'scan_mode'; ... 
'.    Scan resolution  (x, y)            :'    'int   vector'    'scan_resolution'; ... 
'.    Scan percentage                    :'    'int   scalar'    'scan_percentage'; ... 
'.    Recon resolution (x, y)            :'    'int   vector'    'recon_resolution'; ... 
'.    Number of averages                 :'    'int   scalar'    'num_averages'; ... 
'.    Repetition time [msec]             :'    'float scalar'    'repetition_time'; ...   
'.    FOV (ap,fh,rl) [mm]                :'    'float vector'    'fov'; ... 
'.    Slice thickness [mm]               :'    'float scalar'    'slice_thickness'; ...
'.    Slice gap [mm]                     :'    'float scalar'    'slice_gap'; ... 
'.    Water Fat shift [pixels]           :'    'float scalar'    'water_fat_shift'; ... 
'.    Angulation midslice(ap,fh,rl)[degr]:'    'float vector'    'angulation'; ...
'.    Off Centre midslice(ap,fh,rl) [mm] :'    'float vector'    'offcenter'; ... 
'.    Flow compensation <0=no 1=yes> ?   :'    'int   scalar'    'flowcomp'; ...
'.    Presaturation     <0=no 1=yes> ?   :'    'int   scalar'    'presaturation';... 
'.    Cardiac frequency                  :'    'int   scalar'    'card_frequency'; ...
'.    Min. RR interval                   :'    'int   scalar'    'min_rr_interval'; ...
'.    Max. RR interval                   :'    'int   scalar'    'max_rr_interval'; ...
'.    Phase encoding velocity [cm/sec]   :'    'float vector'    'venc'; ... 
'.    MTC               <0=no 1=yes> ?   :'    'int   scalar'    'mtc'; ...
'.    SPIR              <0=no 1=yes> ?   :'    'int   scalar'    'spir'; ...
'.    EPI factor        <0,1=no EPI>     :'    'int   scalar'    'epi_factor'; ...
'.    TURBO factor      <0=no turbo>     :'    'int   scalar'    'turbo_factor'; ...
'.    Dynamic scan      <0=no 1=yes> ?   :'    'int   scalar'    'dynamic_scan'; ...
'.    Diffusion         <0=no 1=yes> ?   :'    'int   scalar'    'diffusion'; ...
'.    Diffusion echo time [msec]         :'    'float scalar'    'diffusion_echo_time'; ...
'.    Inversion delay [msec]             :'    'float scalar'    'inversion_delay'; ...
};