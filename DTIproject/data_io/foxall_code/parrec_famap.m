%
% PARREC_FAMAP
%
% Extract FA maps from PARREC files 
%
%i $Id: parrec_famap.m,v 1.3 2007/01/08 20:59:57 dfoxall Exp dfoxall $
%

filename      = 'H1P31-ARI-HEAD_2_1';
plotname      = 'P31 H1 head coil';
threshold     = 1.0;

%-------------------------------------------------------------------
% READ THE .PAR FILE
% READ THE .REC FILE
% RESCALE & FLOAT THE DATA
%-------------------------------------------------------------------

parfile_name  = [filename,'.PAR'];
recfile_name  = [filename,'.REC'];

[parfile,error_flag]        = read_parfile(parfile_name,'V4');
[int_data,TNI,NLINE,stride] = read_recfile(recfile_name,parfile);
rs_data                     = rescale_parrec(int_data,parfile);

%-------------------------------------------------------------------
% Extract FA Images
%-------------------------------------------------------------------
Nslices                     = parfile.max_slices;
Nsets                       = round(TNI/Nslices);
slice_tag                   = 1;
type_tag                    = 6;
flip_tag                    = 36;

if (Nsets >= 3)
   for I = 1:TNI
       tags  = parfile.scan_tags(I,:);
       slice = tags(slice_tag);
       type  = tags(type_tag);
       flip  = tags(flip_tag); 
  
       if (type == 0)
           I1(slice,:,:) = rs_data(I,:,:); 
       end

       if (type == 2)
           I2(slice,:,:) = rs_data(I,:,:); 
       end

       if (type == 5)
           I3(slice,:,:) = rs_data(I,:,:); 
       end
   end
end

mask  = binary_mask(I2,threshold);
FA    = mask.*I3*flip/100.0;

eval(['save ',filename,' plotname flip parfile I1 I2 I3 FA mask']);
