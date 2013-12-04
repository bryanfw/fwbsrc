


filename = 'IEPI_-_Vanderbilt_WIP_dyn_2_shot_stab_seq_Experiment3_5_1';
plotname = strrep(filename,'_','-');

%----------------------------------------------------------------
% Read Data From PAR/REC files
%----------------------------------------------------------------
parfile_name  = [filename,'.PAR'];
recfile_name  = [filename,'.REC'];

fprintf(1,'%s\n',filename);
[parfile,error_flag]        = read_parfile(parfile_name,'V4');
[int_data,TNI,NLINE,stride] = read_recfile(recfile_name,parfile); 
rs_data                     = rescale_parrec(int_data,parfile);