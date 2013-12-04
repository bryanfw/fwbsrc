%
% READ_PARFILE
%
% Read the scan parameters and scan tag information from 
% a Philips NV .PAR file and return a MATLAB structure
% with the file data stored in named fields.
%
% FORMAT:
%
% [parfile,error_flag] = read_parfile(file_name,version_spec) 
%
% INPUT:
%  
% file_name      string    File path name of the .PAR file 
%
% version_spec   string    M-FILE that defines a version
%                          specification for the general 
%                          information lines contained in  
%                          the file. Supply the M_FILE name 
%                          without the .m extension.
% 
% OUTPUT: 
%
% parfile       structure  A MATLAB structure that contains 
%                          values read from the file in  
%                          named fields according to the M-file
%                          version specification.
%
% error_flag    int        Values > 0 indicate an error was 
%                          encountered trying to fill some
%                          fields inside the parfile defined in 
%                          the version specification. The data 
%                          in the parfile structure may be usable
%                          but is incomplete. 
%
%
% VERSION SPECIFICATION
%
% Philips .PAR files are versioned. The version number of the 
% target file is contained in a comment line such as:
%
% # CLINICAL TRYOUT      Research image export tool  V3 
%
% Indicating the version of the image export tool that wrote
% the target file. Different versions differ in the number of
% general information lines they contain and the length and 
% content of the tag lines that describe the binary image data 
% in the corresponding .REC file. To make this M-FILE more
% flexible, format differences between different versions are
% encoded as a MATLAB cell array in separate M-FILES. The
% experienced user may create his own M-FILE specifications 
% for new versions or subsets of the .PAR file data.
%
% The cell array has the format:
%
% geninfo_spec = { ...
%   line_key_string   value_type_string    field_name_string; ...
%   line_key_string   value_type_string    field_name_string; ...
%   line_key_string   value_type_string    field_name_string; ...
%   line_key_string   value_type_string    field_name_string; ...
% };
%
% Where the line_key_string matches a general information line tag 
% such as:
%
% '.    Examination name                   :'     
%
% and the value_type_string describes the type of data that follows 
% the colon delimeter. Valid value_type_strings are:
%
%  VALUE TYPE            FIELD FORMAT
% --------------------------------------------------------------------- 
% 'char  vector'         format string data as a cell array of white
%                        space separated words.   
% 'char  scaler'         format string data as a single string.    
% 'float vector'         format a numeric data list as an array.
% 'float scaler'         format numeric data as a single value.   
%
% The value_type_strings determine how the information on the rest of
% the specified general information line is formated in the parfile
% structure.
%
% The field_name_strings can be any suitable memnonic for the extracted
% data. The parfile will contain fields whose names are specified by 
% these strings.
%
%
% $Id: read_parfile.m,v 1.3 2006/04/04 15:27:15 dfoxall Exp $
%

function [parfile,error_flag] = read_parfile(file_name,version_spec) 

%-----------------------------------------------------
% INITIALIZATION
%-----------------------------------------------------
error_flag = 0;


%-----------------------------------------------------
% READ IN .PAR A LINE AT A TIME
%
% Establish the number of lines to parse and the
% contents of the file as a cell array of strings 
%
%-----------------------------------------------------
NL  = 0;
fid = fopen(file_name,'r');

if (length(fid) < 1)
    error(['.PAR file ', file_name, ' not found.']);
end

while 1
      curline = fgetl(fid);
      if ( ~ischar(curline) )
          break;
      else
          NL        = NL + 1;
          lines(NL) = cellstr(curline);
      end
end

fclose(fid);


%-----------------------------------------------------
% EXTRACT GENERAL INFORMATION LINES
%
% Establish the number of general information lines
% and copy then into a cell array.
%
%-----------------------------------------------------
NG  = 0;

for L = 1:NL
    curline = char(lines(L));
    cursize = min(size(curline)); 

    if (cursize > 0)
    if (curline(1) == '.')
       NG          = NG + 1;
       geninfo(NG) = cellstr(curline); 
    end
    end

end

if (NG < 1)
   error('.PAR file has invalid format')
end

%-----------------------------------------------------
% EXECUTE THE VERSION SPEC M-FILE
%-----------------------------------------------------
if (length(fopen(version_spec,'r')) > 0)
   eval(version_spec)
   parfile.version = version_spec;
end

if (length(geninfo_spec) < 1)
   error('Invalid version specification');
end


%-----------------------------------------------------
% PARSE GENERAL INFORMATION
%-----------------------------------------------------
NS           = max(size(geninfo_spec)); 
SQ           = sprintf('%c',39);
value_start  = 1 + max(size(char(geninfo_spec(1,1))));

for S=1:NS
    line_key   = char(geninfo_spec(S,1));
    value_type = char(geninfo_spec(S,2));
    field_name = char(geninfo_spec(S,3));
    L          = strmatch(line_key,geninfo);

    if (length(strmatch(line_key,geninfo)) > 0)
        L          = strmatch(line_key,geninfo);
        curline    = char(geninfo(L));
        value_end  = max(size(curline));
    else
        value_type = ':-( VALUE NOT FOUND )-:';
        error_flag = error_flag + 1;
        warning(['Missing line key match. The parfile.',field_name,' field has no useful value.']);
    end 

    switch (value_type)

    case { 'float scalar' 'int   scalar'}
         eval(['parfile.',field_name,' = ',curline(value_start:value_end),';']);

    case { 'float vector' 'int   vector'}
         eval(['parfile.',field_name,' = [',curline(value_start:value_end),'];']);

    case { 'char  scalar' }
         valstr1 = deblank(curline(value_start:value_end));
         valstr2 = strjust(valstr1,'left');
         valstr3 = deblank(valstr2);
         eval(['parfile.',field_name,' = ',SQ,valstr3,SQ,';']); 

    case { 'char  vector' }
         remainder = curline(value_start:value_end);
         ntokens   = 0;
         while (1)
             [token,remainder] = strtok(remainder);
             if ( ~ischar(token) )
                break;
             else
                ntokens         = ntokens + 1;
                tokens(ntokens) = cellstr([' ',SQ,token,SQ,' ']);
             end
             if ( isempty(remainder) )
                break;
             end
         end

         valstr1 = char(tokens(1));
         for T=2:ntokens
             valstr2 = strcat(valstr1,char(tokens(T)));
             valstr1 = valstr2;
         end
         
         eval(['parfile.',field_name,' = {',valstr1,'};']); 

    otherwise
         eval(['parfile.',field_name,' = ',SQ,value_type,SQ,';']); 
    end

end


%-----------------------------------------------------
% EXTRACT SCAN INFORMATION LINES
%-----------------------------------------------------
NI  = 0;

for L = 1:NL
    curline = char(lines(L));
    cursize = min(size(curline)); 

    if (cursize > 0)
    if (curline(1) ~= '.') 
    if (curline(1) ~= '#') 
    if (curline(1) ~= '*') 
       NI           = NI + 1;
       scantags(NI) = cellstr(curline); 
    end
    end
    end
    end

end

if (NI < 1)
   error('Missing scan information in .PAR file');
end

for I=1:NI
    eval(['parfile.scan_tags(',num2str(I),',:) = [',char(scantags(I)),'];']);
end
