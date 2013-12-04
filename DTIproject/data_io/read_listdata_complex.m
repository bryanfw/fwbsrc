function [complex_img] = read_listdata_complex(lstfile)
% Read the complex data file obtained from a Philips scanner. The file name
% is something like cpx_XXX.list, cpx_XXX.data
%
% Inputs:
%           If no inputs, via UI, select name of .list file
%           assumes that the corresponding .data file is in the same folder
%           with .list file.
%
% Output: complex_img
%
%==========================================================================
% Author        Hui Wang, Ph.D student
%               University of Louisville, KY
%               12/20/2010
% Updated: 8/15/2011
%==========================================================================
tic;
if nargin < 1  % Select file from UI
    [fname,pname] = uigetfile('*.list','Select *.list file');
    if isequal(fname, 0)
        disp('User selection has been cancelled.');
    else
        lstfile=[pname fname];
        [lstparms, DVattribs] = read_listfile_complex(lstfile);
        [complex_img] = read_datafile_complex([lstfile(1:(end-5)) '.data'], lstparms, DVattribs);
    end
else          % Select file from input arguments
    [lstparms, DVattribs] = read_listfile_complex(lstfile);
    [complex_img] = read_datafile_complex([lstfile(1:(end-5)) '.data'], lstparms, DVattribs);
end
toc;

end
