function print_status_bar(idx, num, varargin)
% print_status_bar - prints a simple status bar for updating the user
%
% Three forms:
% 1) print_status_bar(i, num)
% 2) print_status_bar(i, num, num_brks)
% 3) print_status_bar(i, num, num_brks, pre_str)
%
% Input: idx - the current index into the full list
%        num - the total length of the list
%        num_brks - the number of status updates to provide
%        pre_str - the constant string appearing out front
%
% Output: (None)

if length(varargin) == 0
    if (num < 78)
        num_brks = num;
    else
        num_brks = 78;
    end
    pre_str = '';
elseif length(varargin) == 1
    if isstr(varargin{1})
        pre_str = varargin{1};
        if (num < 78)
            num_brks = num;
        else
            num_brks = 78;
        end
    else
        num_brks = varargin{1};
        pre_str = '';
    end
elseif length(varargin) == 2
    num_brks = varargin{1};
    pre_str = varargin{2};
else
    error('Too many input arguments');
end

cval = floor((num_brks * idx-1) / num)+1;
pval = floor((num_brks * (idx-2)) / num)+1;

rmstr = '';
if (cval > pval)
    str = ['[', repmat('=', [1 cval]), repmat('_', [1 num_brks-cval]), ']'];
    if (idx > 1)
        rmstr = repmat('\b', [1 num_brks+2]);
    elseif (idx == 1)
        rmstr = pre_str;
    end
    fprintf([rmstr, str]);
end

if (idx == num)
    fprintf('\n');
end
