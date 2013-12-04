%
% READ DATA
%
% Read the data from a .data file obtained on a Philips scanner. In this version PHX and FRX
% and PHC and FRC are read into the data structure if they are encountered.
%
% SYNOPSIS:
%
% sort_data = read_data(listfile,data_filename);
%
% INPUT:
%
% listfile      (structure)    MATLAB structure extracted from
%                              .LIST file by read_list.m
%
% data_filename (string)       Name of .data file with .data
%                              extension.
%
% OUTPUT
%
% sort_data     (structure)    The data lines extracted from
%                              the .data file sorted by type
%                              along with their sizes, counts
%                              and data description attributes.
%
% sort_data =
%          noi: [1x40000 double]
%      noi_att: [0 0 0 0 0 26 0 0 0 0 0 0 1 0 0 0 0 0 160000 0]
%          std: [4800x800 double]
%      std_att: [4800x20 double]
%    std_count: 4800
%     std_size: 800
%    noi_count: 1
%     noi_size: 40000
%    nav_count: 0
%     nav_size: 0
%    frc_count: 0
%     frc_size: 0
%    phc_count: 0
%     phc_size: 0
%    rej_count: 0
%     rej_size: 0
%
% $Id:$
%

function sort_data = read_data(listfile,data_filename);

%--------------------------------------------------------------------
% CONSTANTS
%
%--------------------------------------------------------------------
FLOAT32_BYTES = 4;


%--------------------------------------------------------------------
% INITIALIZATION
%
% dtc                  (structure)   Data type count
% ddi                  (structure)   Data description indices
%
%--------------------------------------------------------------------
error_flag  = 0;
dtc.std_cnt = 0;
dtc.std_stp = 0;
dtc.rej_cnt = 0;
dtc.rej_stp = 0;
dtc.phx_cnt = 0;
dtc.phx_stp = 0;
dtc.phc_cnt = 0;
dtc.phc_stp = 0;
dtc.frx_cnt = 0;
dtc.frx_stp = 0;
dtc.frc_cnt = 0;
dtc.frc_stp = 0;
dtc.noi_cnt = 0;
dtc.noi_stp = 0;
dtc.nav_cnt = 0;
dtc.nav_stp = 0;

data_attribute_size = size(listfile.data_attributes);
num_attributes      = data_attribute_size(2);

switch (num_attributes)
       case 19
            ddi_19keys
       case 20 
            ddi_20keys
       otherwise
            error('Check number of data attributes in .LIST file'); 
end


%--------------------------------------------------------------------
% CALCULATE DATA FILE SIZE FROM THE DATA DESCRIPTION
% READ THE .data FILE (IEEE FLOATS) 
%--------------------------------------------------------------------
attribute_size   = size(listfile.data_attributes);
total_calc_bytes = 0;

for A = 1:attribute_size(1)
    total_calc_bytes = total_calc_bytes + listfile.data_attributes(A,ddi.size); 
end

total_calc_floats        = total_calc_bytes/FLOAT32_BYTES;
fid                      = fopen(data_filename,'r');
[raw_data,raw_data_size] = fread(fid,inf,'float32'); 
fclose(fid);


%--------------------------------------------------------------------
% COUNT DATA LINES BY LINE TYPE 
% COPY SORT DATA LINES BY TYPE INTO OUTPUT STRUCTURE 
%--------------------------------------------------------------------
copy_beg = 1;
copy_end = 0;
copy_off = 0;

for A = 1:attribute_size(1)
    copy_beg = copy_beg + copy_off;
    copy_off = listfile.data_attributes(A,ddi.size)/FLOAT32_BYTES;
    copy_end = (copy_beg - 1) + copy_off;
    copy_typ = char(listfile.data_type(A));

    switch ( copy_typ )
    case { 'NOI' }
            dtc.noi_cnt = dtc.noi_cnt + 1;
            if (dtc.noi_cnt == 1)
                dtc.noi_stp = listfile.data_attributes(A,ddi.size)/FLOAT32_BYTES;
            end
            sort_data.noi(dtc.noi_cnt,:)     = raw_data(copy_beg:copy_end);
            sort_data.noi_att(dtc.noi_cnt,:) = listfile.data_attributes(A,:);

    case { 'STD' }
            dtc.std_cnt = dtc.std_cnt + 1;
            if (dtc.std_cnt == 1)
                dtc.std_stp = listfile.data_attributes(A,ddi.size)/FLOAT32_BYTES;
            end
            sort_data.std(dtc.std_cnt,:)     = raw_data(copy_beg:copy_end);
            sort_data.std_att(dtc.std_cnt,:) = listfile.data_attributes(A,:);

    case { 'REJ' }
            dtc.rej_cnt = dtc.rej_cnt + 1;
            if (dtc.rej_cnt == 1)
                dtc.rej_stp = listfile.data_attributes(A,ddi.size)/FLOAT32_BYTES;
            end
            sort_data.rej(dtc.rej_cnt,:)     = raw_data(copy_beg:copy_end);
            sort_data.rej_att(dtc.rej_cnt,:) = listfile.data_attributes(A,:);

    case { 'PHX' }
            dtc.phx_cnt = dtc.phx_cnt + 1;
            if (dtc.phx_cnt == 1)
                dtc.phx_stp = listfile.data_attributes(A,ddi.size)/FLOAT32_BYTES;
            end
            sort_data.phx(dtc.phx_cnt,:)     = raw_data(copy_beg:copy_end);
            sort_data.phx_att(dtc.phx_cnt,:) = listfile.data_attributes(A,:);

    case { 'PHC' }
            dtc.phc_cnt = dtc.phc_cnt + 1;
            if (dtc.phc_cnt == 1)
                dtc.phc_stp = listfile.data_attributes(A,ddi.size)/FLOAT32_BYTES;
            end
            sort_data.phc(dtc.phc_cnt,:)     = raw_data(copy_beg:copy_end);
            sort_data.phc_att(dtc.phc_cnt,:) = listfile.data_attributes(A,:);

    case { 'FRX' }
            dtc.frx_cnt = dtc.frx_cnt + 1;
            if (dtc.frx_cnt == 1)
                dtc.frx_stp = listfile.data_attributes(A,ddi.size)/FLOAT32_BYTES;
            end
            sort_data.frx(dtc.frx_cnt,:)     = raw_data(copy_beg:copy_end);
            sort_data.frx_att(dtc.frx_cnt,:) = listfile.data_attributes(A,:);

    case { 'FRC' }
            dtc.frc_cnt = dtc.frc_cnt + 1;
            if (dtc.frc_cnt == 1)
                dtc.frc_stp = listfile.data_attributes(A,ddi.size)/FLOAT32_BYTES;
            end
            sort_data.frc(dtc.frc_cnt,:)     = raw_data(copy_beg:copy_end);
            sort_data.frc_att(dtc.frc_cnt,:) = listfile.data_attributes(A,:);

    case { 'NAV' }
            dtc.nav_cnt = dtc.nav_cnt + 1;
            if (dtc.nav_cnt == 1)
                dtc.nav_stp = listfile.data_attributes(A,ddi.size)/FLOAT32_BYTES;
            end
            sort_data.nav(dtc.nav_cnt,:)     = raw_data(copy_beg:copy_end);
            sort_data.nav_att(dtc.nav_cnt,:) = listfile.data_attributes(A,:);

    otherwise
         error(['read_data: Unexpected data line type encountered: <',copy_typ,'>']);
    end 

end
 
sort_data.std_count = dtc.std_cnt;
sort_data.std_size  = dtc.std_stp;
sort_data.noi_count = dtc.noi_cnt;
sort_data.noi_size  = dtc.noi_stp;
sort_data.nav_count = dtc.nav_cnt;
sort_data.nav_size  = dtc.nav_stp;
sort_data.frx_count = dtc.frx_cnt;
sort_data.frx_size  = dtc.frx_stp;
sort_data.frc_count = dtc.frc_cnt;
sort_data.frc_size  = dtc.frc_stp;
sort_data.phx_count = dtc.phx_cnt;
sort_data.phx_size  = dtc.phx_stp;
sort_data.phc_count = dtc.phc_cnt;
sort_data.phc_size  = dtc.phc_stp;
sort_data.rej_count = dtc.rej_cnt;
sort_data.rej_size  = dtc.rej_stp;

%-------------
% END
%-------------


