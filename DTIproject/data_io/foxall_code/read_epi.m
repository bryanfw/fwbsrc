%
% READ_EPI - Read in an raw epi data set for processing using the read_list/read_data macros 
%          
%

%---------------------------------------------------------------
% Read the data description in the .LIST file
%---------------------------------------------------------------
filename              = 'raw_600';
plotname              = strrep(filename,'_','-');
listfile              = [filename,'.list'];
datafile              = [filename,'.data'];
[data_des,error_flag] = read_list(listfile,'SSEPI');


%---------------------------------------------------------------
% Read the raw data from the .DATA file
%---------------------------------------------------------------
data_raw              = read_data(data_des,datafile);


