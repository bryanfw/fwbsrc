%
% READ_LIST
%
% Read the data description parameters in an exported .list file 
% obtained from a Philips scanner.
%
% FORMAT:
%
% [listfile,error_flag] = read_list(file_name,conversion_spec);
%
%
% INPUT:
%
% file_name           (string)      File path name of the .list file
% 
% conversion_spec     (string)      M-FILE that defines the conversion
%                                   specification for the .list information
%                                   Supply the M-FILE name without the .m
%                                   extension.
%
% OUTPUT:
%
% listfile            (structure)   A MATLAB structure that contains data
%                                   description information extracted from
%                                   the .list file according to the conversion
%                                   specification.
%
% error_flag          int           Values greater than 0 indicate an error
%                                   was encountered tyring to extract
%                                   data description information and fill                                
%                                   the fields expected in the conversion
%                                   specification.
%
%
% CONVERSION SPECIFICATION
%
% An M-FILE is used to provide a conversion specifcation for the data 
% description information. This allows the basic .list file reader
% to be customized and adapted for different applications. Conversion  
% specifications are provided for both general information lines
% and for comment lines in the form of MATLAB cell arrays
%
% The cell arrays have the format:
%
% comment_info_spec = { ...
%   line_key_string   value_type_string    field_name_string; ...
%   line_key_string   value_type_string    field_name_string; ...
% };
%
% geninfo_spec = { ...
%   line_key_string   value_type_string    field_name_string; ...
%   line_key_string   value_type_string    field_name_string; ...
% };
%
% Where the line_key_string matches a comment line or
% general data information line tag such as:
%
% '#    Scan name :'     
%
% OR
%
% '.    0  0  0  number_of_mixes                  :'     
%
% The value_type_string describes the type of data that follows 
% the colon delimeter. Valid value_type_strings are:
%
%  VALUE TYPE            FIELD FORMAT
% --------------------------------------------------------------------- 
% 'char  vector'         format string data as a cell array of white
%                        space separated words.   
% 'char  scaler'         format string data as a single string.    
% 'float vector'         format a floating data list as an array.
% 'float scaler'         format a single floating value.   
% 'int   vector'         format an integer data list as an array.
% 'int   scaler'         format a single integer value.   
%
% The value_type_strings determine how the information on the rest of
% the specified general information line is formated in the parfile
% structure.
%
% The field_name_strings can be any suitable memnonic for the extracted
% data. The parfile will contain fields whose names are specified by 
% these strings.
%
% $Id: read_list.m,v 1.2 2006/05/25 14:24:43 dfoxall Exp $
%

function [listfile,error_flag] = read_list(file_name,conversion_spec);

%-----------------------------------------------------
% INITIALIZATION
% EXECUTE THE VERSION SPEC M-FILE
%-----------------------------------------------------
error_flag = 0;

if (length(fopen(conversion_spec,'r')) > 0)
   eval(conversion_spec)
   listfile.conversion_spec = conversion_spec;
end

if (length(comment_info_spec) < 1)
   error('Invalid conversion specification');
end

if (length(geninfo_spec) < 1)
   error('Invalid conversion specification');
end



%-----------------------------------------------------
% READ IN THE FILE A LINE AT A TIME
%
% Establish the number of lines to parse and the
% contents of the file as a cell array of strings 
%
%-----------------------------------------------------
NL  = 0;
fid = fopen(file_name,'r');

if (length(fid) < 1)
    error(['.list file ', file_name, ' not found.']);
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
% EXTRACT COMMENT LINES
%
% Establish the number of comment lines
% and copy then into a cell array.
%
%-----------------------------------------------------
NC  = 0;

for L = 1:NL
    curline = char(lines(L));
    cursize = min(size(curline)); 

    if (cursize > 0)
    if (curline(1) == '#')
       NC               = NC + 1;
       comment_info(NC) = cellstr(curline); 
    end
    end

end

if (NC < 1)
   error('.list file has invalid format')
end



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
   error('.list file has invalid format')
end




%-----------------------------------------------------
% PARSE COMMENT INFORMATION
%-----------------------------------------------------
NS           = max(size(comment_info_spec)); 
SQ           = sprintf('%c',39);

for S=1:NS
    line_key     = char(comment_info_spec(S,1));
    value_type   = char(comment_info_spec(S,2));
    field_name   = char(comment_info_spec(S,3));
    L            = strmatch(line_key,comment_info);
    value_start  = 1 + max(size(line_key));

    if (length(strmatch(line_key,comment_info)) > 0)
        L          = strmatch(line_key,comment_info);
        curline    = char(comment_info(L));
        value_end  = max(size(curline));
    else
        value_type = ':-( VALUE NOT FOUND )-:';
        error_flag = error_flag + 1;
        warning(['Missing line key match. The listfile.',field_name,' field has no useful value.']);
    end 

    switch (value_type)

    case { 'float scalar' 'int   scalar'}
         eval(['listfile.',field_name,' = ',curline(value_start:value_end),';']);

    case { 'float vector' 'int   vector'}
         eval(['listfile.',field_name,' = [',curline(value_start:value_end),'];']);

    case { 'char  scalar' }
         valstr1 = deblank(curline(value_start:value_end));
         valstr2 = strjust(valstr1,'left');
         valstr3 = deblank(valstr2);
         eval(['listfile.',field_name,' = ',SQ,valstr3,SQ,';']); 

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
             if ( ~ischar(remainder) )
                break;
             end
         end

         valstr1 = char(tokens(1));
         for T=2:ntokens
             valstr2 = strcat(valstr1,char(tokens(T)));
             valstr1 = valstr2;
         end
         
         eval(['listfile.',field_name,' = {',valstr1,'};']); 

    otherwise
         eval(['listfile.',field_name,' = ',SQ,value_type,SQ,';']); 
    end

end



%-----------------------------------------------------
% PARSE GENERAL INFORMATION
%-----------------------------------------------------
NS           = max(size(geninfo_spec)); 
SQ           = sprintf('%c',39);
%value_start  = 1 + max(size(char(geninfo_spec(1,1))));

for S=1:NS
    line_key     = char(geninfo_spec(S,1));
    value_type   = char(geninfo_spec(S,2));
    field_name   = char(geninfo_spec(S,3));
    L            = strmatch(line_key,geninfo);
    value_start  = 1 + max(size(line_key));

    if (length(strmatch(line_key,geninfo)) > 0)
        L          = strmatch(line_key,geninfo);
        curline    = char(geninfo(L));
        value_end  = max(size(curline));
    else
        value_type = ':-( VALUE NOT FOUND )-:';
        error_flag = error_flag + 1;
        warning(['Missing line key match. The listfile.',field_name,' field has no useful value.']);
    end 

    switch (value_type)

    case { 'float scalar' 'int   scalar'}
         eval(['listfile.',field_name,' = ',curline(value_start:value_end),';']);

    case { 'float vector' 'int   vector'}
         eval(['listfile.',field_name,' = [',curline(value_start:value_end),'];']);

    case { 'char  scalar' }
         valstr1 = deblank(curline(value_start:value_end));
         valstr2 = strjust(valstr1,'left');
         valstr3 = deblank(valstr2);
         eval(['listfile.',field_name,' = ',SQ,valstr3,SQ,';']); 

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
             if ( ~ischar(remainder) )
                break;
             end
         end

         valstr1 = char(tokens(1));
         for T=2:ntokens
             valstr2 = strcat(valstr1,char(tokens(T)));
             valstr1 = valstr2;
         end
         
         eval(['listfile.',field_name,' = {',valstr1,'};']); 

    otherwise
         eval(['listfile.',field_name,' = ',SQ,value_type,SQ,';']); 
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
   error('Missing scan information in .list file');
end

for I=1:NI
    curline = char(scantags(I));
    [token,remainder] = strtok(curline);
    eval(['listfile.data_attributes(',num2str(I),',:) = [',char(remainder),'];']);
    eval(['listfile.data_type(',num2str(I),') = cellstr(token);']);
end
