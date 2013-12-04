
function [s] = f_cleanFieldname(s)
%[f_cleanFieldname] cleans unnecessary characters reading field name
%in .LIST file.
%
%
% Last modified
% 2010.12.16
%   Generated.



%% Main
illegal_chars = {...
    '+','-','*','.',...
    '^','\','/','.',...
    '=','~','<','>',...
    '&','|',':',';',...
    '(',')','{','}',...
    '[',']','{','}',...
    '''','%',' ','!', ...
    '@','#','$','`',...
    '?',',','"',...
    };

general_replacement_char = '_';
firstchar_replacement_char = 'x'; % cannot be an underscore

for k=1:length(illegal_chars),
    s = strrep(s,illegal_chars{k},general_replacement_char);
end

% first character cannot be a number
firstchar_code = double(s(1));
if ( (firstchar_code>=double('0')) && (firstchar_code<=double('9')) )
    s(1) = firstchar_replacement_char;
end

% first character cannot be an underscore
if(s(1)=='_'),
    s(1) = firstchar_replacement_char;
end



%% END






