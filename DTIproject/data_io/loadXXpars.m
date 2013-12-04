%% LOADXXPARS     Return a structure of information contained within Philips DIXOM XX_ file
%
% [XXPARS] = LOADXXPARS(XX_FILENAME)
%
%   XXPARS is a structure containing parameters stored within certain of
%   the XX_ files accompanying Philips DICOM data sets. Many useful
%   parameter values are contained within such XX_ files
%
%  See also: LOADDICOM
%
%  Dependencies: DICOM-DICT-PHILIPS.TXT
%

%% Revision History
% * 2009.03.20    initial version - brianwelch

function XXpars = loadXXpars(xx_filename)

%% Open XX file and read "pixel data"
dicomdict('set','dicom-dict-philips.txt');
info = dicominfo(xx_filename);

XXpars = [];
itemnames = fieldnames(info.MRBlobDataObjectArray);
for k=1:length(itemnames),
    MRBlob_name = info.MRBlobDataObjectArray.(itemnames{k}).MRBlobName;
    XXpars.(MRBlob_name).MRTypeName = info.MRBlobDataObjectArray.(itemnames{k}).MRTypeName;
    XXpars.(MRBlob_name).parameters = parseMRBlobData(info.MRBlobDataObjectArray.(itemnames{k}).MRBlobData);
end

function parameter_struct = parseMRBlobData(MRBlobData)
tmpfile = sprintf('%s/parseMRBlobData.hex', ctfroot);
fid = fopen(tmpfile,'wb');
fwrite(fid,MRBlobData,'int16');
fclose(fid);

protocol_header = protocol_header_info(tmpfile);
parameter_struct = protocol_parameters(tmpfile, protocol_header.number_of_parameters);
parameter_struct = get_protocol_parameters_values(parameter_struct,tmpfile);

% #
% # Subroutine for unpacking protocol_buffer header
% #  
% # /* Structure from Philips source code */
% #typedef struct
% #{
% #    SRP		int_pool_srp;		/* int []		*/
% #    USHORT		number_of_integers;
% #    USHORT		struct_type;		/* Type of this struct	*/
% #						/* GGPARS_STRUCT_TYPE_ENUM*/
% #						/* Must be 0		*/
% #						/* MAY NOT BE MOVED!	*/
% #    SRP		float_pool_srp;		/* float []		*/
% #    UINT		number_of_floats;
% #    SRP		string_pool_srp;	/* GGPARS_STRING_ADT []	*/
% #    UINT		number_of_bytes;
% #    SRP		parameters_srp;		/* GGPARS_PARAMETER_STRUCT []*/
% #    UINT		number_of_parameters;
% #} GGPARS_INDEX_STRUCT;
% #
function info = protocol_header_info(tmpfile),

unpackstr = 'iSSiIiIiI';
unpacknames = {'int_pool_srp','number_of_integers','struct_type','float_pool_srp','number_of_floats','string_pool_srp','number_of_bytes','parameters_srp','number_of_parameters'};

info = [];
fid = fopen(tmpfile,'r');
for k=1:length(unpackstr),
    info.(unpacknames{k}) = fread(fid,readlength(unpackstr,k),readtype(unpackstr,k)); 
end
fclose(fid);

% #
% # Subroutine for unpacking protocol_buffer parameters
% #  
% # /* Structure from Philips source code */
% #typedef struct
% #{
% #    char			par_name[ GGPARS_MAX_ID_STRING_1 ]; /* 32 + 1 */
% #    BOOLEAN_8		par_enab_for_disp;
% #    int			par_type;
% #    UINT			par_dimension;
% #    UINT			par_level;
% #    SRP			par_value_srp;	/* GGPARS_VALUE_UNION		*/
% #} GGPARS_PARAMETER_STRUCT;
% #
function parameters = protocol_parameters(tmpfile, number_of_parameters)

headerunpackstr = 'iSSiIiIiI';
headerlength = readlengthtotal(headerunpackstr);

unpackstr = 'A33ciIIi';
unpacknames = {'par_name','par_enab_for_disp','par_type','par_dimension','par_level','par_value_srp'};
recordlength = readlengthtotal(unpackstr);

parameters=[];
fid = fopen(tmpfile,'r');
fseek(fid,headerlength,'bof');
for n=1:number_of_parameters,    
    for k=1:length(unpacknames),
        if k==1,
            tmp = fread(fid,[1 readlength(unpackstr,k)],readtype(unpackstr,k));
            toks = regexp(tmp,'([a-z_A-Z0-9]+)','tokens');
            par_name = char(toks{1});
        else
            parameters.(par_name).(unpacknames{k}) = fread(fid,[1 readlength(unpackstr,k)],readtype(unpackstr,k));
        end
    end
end
fclose(fid);

% # 
% #	typedef enum
% #	{
% #		GGPARS_VALUE_T_MIN = -1,
% #		GGPARS_VALUE_T_FLOAT,
% #		GGPARS_VALUE_T_INTEGER,
% #		GGPARS_VALUE_T_STRING,
% #		GGPARS_VALUE_T_TGM,		/* obsolete; to be removed */
% #		GGPARS_VALUE_T_ENUM,
% #		GGPARS_VALUE_T_MAX
% #	} GGPARS_VALUE_TYPE_ENUM;
% #
function parameter_struct = get_protocol_parameters_values(parameter_struct,tmpfile)

fid = fopen(tmpfile,'r','ieee-le');
%fid = fopen(tmpfile,'r','ieee-be');

headerlength = readlengthtotal('iSSiIiIiI');
recordlength = readlengthtotal('A33ciIIi');

par_names = fieldnames(parameter_struct);
for k=1:length(par_names),
    par_type = parameter_struct.(par_names{k}).par_type;
    par_dimension = parameter_struct.(par_names{k}).par_dimension;
    par_value_srp = parameter_struct.(par_names{k}).par_value_srp;
    
    valuesPtr = headerlength + k*recordlength + par_value_srp - 4;
    fseek(fid,valuesPtr,'bof');
    
    switch par_type,
        case 1,
            par_type_name = 'int';
            unpackstr = sprintf('i%d',par_dimension);
        case 4,
            par_type_name = 'enum';
            unpackstr = sprintf('i%d',par_dimension);
        case 0,
            par_type_name = 'float';
            unpackstr = sprintf('f%d',par_dimension);
        case 2,
            par_type_name = 'string';
            unpackstr = sprintf('c%d',par_dimension);
        otherwise
            error(sprintf('Unknown par_type: (%s)',par_type));
    end

    parameter_struct.(par_names{k}).par_type_name = par_type_name;
    parameter_struct.(par_names{k}).value = fread(fid,[1 readlength(unpackstr,1)],readtype(unpackstr,1));
    
    if strcmp(par_type_name,'int');
        idx = find(parameter_struct.(par_names{k}).value>=2147450879);
        parameter_struct.(par_names{k}).value(idx) = 0;
    end
end

function readtype_str = readtype(unpackstr,k),

readinfo = regexp(unpackstr,'(?<readtype_char>[AciISf])(?<readtype_cnt>\d*)','names');

switch char(readinfo(k).readtype_char),
    case 'A',
        readtype_str = 'uint8=>char';
    case 'c',
        readtype_str = 'int8';
    case 'S',
        readtype_str = 'uint16';
    case 'i',
        readtype_str = 'int32';
    case 'I',
        readtype_str = 'uint32';
    case 'f',
        readtype_str = 'float32';
end

function readtype_cnt = readlength(unpackstr,k),

readinfo = regexp(unpackstr,'(?<readtype_str>[AciISf])(?<readtype_cnt>\d*)','names');
readtype_cnt = str2num(readinfo(k).readtype_cnt);
if isempty(readtype_cnt),
    readtype_cnt = 1;
end


function rl = readlengthtotal(unpackstr),

readinfo = regexp(unpackstr,'(?<readtype_str>[AciIS])(?<readtype_cnt>\d*)','names');

rl = 0;
for k=1:length(readinfo),
    readtype_str = char(readinfo(k).readtype_str);
    readtype_cnt = str2num(readinfo(k).readtype_cnt);
    if isempty(readtype_cnt),
        readtype_cnt = 1;
    end
    switch readtype(unpackstr,k),
        case {'int8','uint8=>char'},
            rl = rl + 1*readtype_cnt;
        case {'uint16'},
            rl = rl + 2*readtype_cnt;
        case {'int32','uint32','float32'},
            rl = rl + 4*readtype_cnt;
    end
end