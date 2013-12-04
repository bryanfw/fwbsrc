
function k_corr_m = f_epi_phase_corr_phc(k_m,phc_data,DATAFLAGparams, ...
    DWIparams,ind_sl,ind_ec,ind_coil,ind_avg,ind_dw)
%[f_epi_phase_corr_phc] does EPI echo phase correction using PHC data.
%
% USAGE
%   k_corr_m = f_epi_phase_corr_phc(k_m,phc_data,DATAFLAGparams, ...
%     DWIparams,ind_sl,ind_ec,ind_coil,ind_dw,ind_avg)
%
% INPUT
%   k_m:        k-space STD data for each slice (then it is 2D)
%   phc_data:   PHC data
%   DWIparams:  Following field values are used here,
%       nSHOT (numebr of shot)
%       EPI_FACTOR (EPI factor for image- and navigator-echo)
%       STD_start_sign (starting gradient polarity of STD data,
%           'sign' attribute in .LIST)
%       PHC_start_sign (starting gradient polarity of PHC data,
%           'sign' attribute in .LIST)
%   DATAFLAGparams: Following field values are used here,
%       epiCorrFitMethod (linear phase fitting method used in
%           here)
%   ind_sl:     index for slice
%   ind_ec:     index for echo (image or nagivator)
%   ind_coil:   index for coil
%   ind_avg:    index for average (works for single-shot, b=0, ind_avg=2 case)
%   ind_dw:     index for diffusion weighting or not, [0,1]
%
% OUTPUT
%   k_corr_m:   Corrected k-space data
%
%
% Last modified
% 2010.07.06.
% 2010.07.07. Modified for MSH also.
% 2010.07.13. Add flag_fitPhase == 4 method.
%             flag_fitPhase==1 is the best for msh-EPI echo phase
%             correction.
% 2010.07.26. 'DATAFLAGparams.epiCorrFitMethod' is taking 'flag_fitPhase'
%             as an input.
% 2010.08.09.
%   This is generated from [epi_phase_corr_phc_v3.m].
%   This contains all comments.
% 2010.10.22.
%   Take care of DWIparams.PHC_start_sign.
%   Take 'ind_dw' and 'ind_avg' as an input for correcting single-shot,
%   b=0, ind_avg=2 case, and take 'conj(phase_mult_m)' for this case.



%% Check input

% Set flag for running while statement only once.
%* flag_fitPhase==1: Fit to average of 'pi_subtr_m' along phase-encoding
%                    direction. Default.
%* flag_fitPhase==2: Fit to each ky line of 'pi_subtr_m'.
%* flag_fitPhase==3: Do Not Use This.
%* flag_fitPhase==4: Don't fit. Use the whole 'pi_subtr_m' for phase
%                    correction.
%* Usually flag_fitPhase==1 generates better EPI echo phase correction and
%* recon images than using others.

% if DWIparams.nSHOT == 1
%     flag_fitPhase = 2; % default for single-shot
% else
%     %* flag_fitPhase==1 is better for image-echo than 2.
%     %* flag_fitPhase==1 is same or better for navigator-echo as 2.
%     %* flag_fitPhase==4 is not good as flag_fitPhase==1.
%     flag_fitPhase = 4; % default for multi-shot is 1
% end
flag_fitPhase = DATAFLAGparams.epiCorrFitMethod;
flag_temp = flag_fitPhase;

% Get the number of ky lines, EPI factor * Nshot.
ny = DWIparams.nSHOT * DWIparams.EPI_FACTOR(ind_ec);

% Get Nx.
Nx = size(k_m,2);



%% Calculate the phase difference between +ve and -ve EPI echo phase

% Reserve output.
k_corr_m = k_m * 0;


% Fixed indices in PHC.
idx_ky = 1;
idx_kz = 1;
idx_dyn = 1;
idx_card = 1;
idx_mix = 1;


% Read PHC data in +ve and -ve gradient for each shot.
% p1_m = squeeze(phc_data(:,idx_ky,idx_kz, ...
%     ind_sl,idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,1,:)).';
% p2_m = squeeze(phc_data(:,idx_ky,idx_kz, ...
%     ind_sl,idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,2,:)).';

%* PHC has 2 averages, but .LIST file does not denote 'aver' attribute.
%* 'aver==1 or 2' uses just one average without reading 'aver' attribute
%* and 'aver==3' uses 'aver' and read each 'aver' data. It doesn't make
%* much difference.
aver = 3;
if aver==1 || aver==2    
    p1_m = squeeze(phc_data(:,idx_ky,idx_kz, ...
        ind_sl,idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,aver,1,:)).';
    p2_m = squeeze(phc_data(:,idx_ky,idx_kz, ...
        ind_sl,idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,aver,2,:)).';
elseif aver==3
    p1_m = squeeze(mean(phc_data(:,idx_ky,idx_kz, ...
        ind_sl,idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,:,1,:),10)).';
    p2_m = squeeze(mean(phc_data(:,idx_ky,idx_kz, ...
        ind_sl,idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,:,2,:),10)).';
    
    %----- TEMP:2010.09.17 -----
%     phc11_m = squeeze(phc_data(:,idx_ky,idx_kz, ...
%         ind_sl,idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,1,1,:)).';
%     phc12_m = squeeze(phc_data(:,idx_ky,idx_kz, ...
%         ind_sl,idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,2,1,:)).';
%     phc21_m = squeeze(phc_data(:,idx_ky,idx_kz, ...
%         ind_sl,idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,1,2,:)).';
%     phc22_m = squeeze(phc_data(:,idx_ky,idx_kz, ...
%         ind_sl,idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,2,2,:)).';
%     p1_m = (phc11_m+conj(phc12_m))/2;
%     p2_m = (phc21_m+conj(phc22_m))/2;
    %----- TEMP: 2010.09.17 -----
end


% Generate ky index for +ve and -ve gradient in PHC.
%* Though STD data has, e.g., 4 shots with 17 EPI factors,
%* PHC data has only 17 'grad' property meaning that
%* only 9 +ve and 8 -ve gradient PHC data are enough
%* for EPI echo phase correction.
ind_odd_ky = 1:2:DWIparams.EPI_FACTOR(ind_ec);
ind_even_ky = 2:2:DWIparams.EPI_FACTOR(ind_ec);


% Read PHC data in +ve and -ve gradient.
%* This will make even ky lines are always smaller or equal to odd ky lines.
if DWIparams.PHC_start_sign==1
    pk_odd_m = p1_m(ind_odd_ky,:);
    pk_even_m = p2_m(ind_even_ky,:);
else
    pk_odd_m = p2_m(ind_odd_ky,:);
    pk_even_m = p1_m(ind_even_ky,:);
end


% Get image-space data.
pi_odd_m = ift1(pk_odd_m,2);
pi_even_m = ift1(pk_even_m,2);


% Calculate the length of ky index which is multiplied to even data.
% Only ind_ind_ky_mult = 1:length(ind_ky_even) is enough, but this is
% to show explicitly.
if mod(DWIparams.EPI_FACTOR(ind_ec),2)~=0 % odd has one more ky line
    ind_ind_ky_mult = 1:size(ind_even_ky,2);
else
    ind_ind_ky_mult = 1:size(ind_even_ky,2);
end


% Subtract even phase from odd image-space phase.
pi_subtr_m = pi_odd_m(ind_ind_ky_mult,:) .* conj(pi_even_m);
if strcmpi(DWIparams.shot_mode,'single-shot') && ...
        ind_dw==0 && ind_avg==2
    pi_subtr_m = pi_even_m .* conj(pi_odd_m(ind_ind_ky_mult,:));
end


% Normalize and make NaN and Inf zero.
pi_subtr_m = pi_subtr_m./abs(pi_subtr_m);
pi_subtr_m(isnan(pi_subtr_m)) = 0;
pi_subtr_m(isinf(pi_subtr_m)) = 0;



%==========================================
% Fit to ky-directionally averaged phase.
% flag_fitPhase = 1.
while flag_temp==1
    
    %* Average phase over ky direction.
    pi_subtr_avg_v = sum(pi_subtr_m,1)/size(pi_subtr_m,1);
    %phase_v = pi_subtr_avg_v; % use phase_v after fitting to linear
    
    
    %* Take only central phase values from pi_subtr_avg_v.
    center_length = length(pi_subtr_avg_v)/4;
    take_center = floor(length(pi_subtr_avg_v)/2)+1;
    start_left = take_center - floor(center_length/2);
    end_right = take_center + floor(center_length/2)-1;
    pi_subtr_avg_cut_v = pi_subtr_avg_v*0;
    pi_subtr_avg_cut_v(start_left:end_right) = ...
        pi_subtr_avg_v(start_left:end_right);
    %** Test plot.
    %figure,plot(1:Nx,angle(pi_subtr_avg_v),'-r',...
    %  1:Nx,angle(pi_subtr_avg_cut_v),':b')
    %title(sprintf('pi\\_subtr\\_avg\\_v(r), pi\\_subtr\\_avg\\_cut\\_v(b)'))
    
    
    %* Fit the central part of phase to linear function.
    phase_unwrap_v = unwrap(angle(pi_subtr_avg_cut_v)); % unwrap before fitting
    p = polyfit(start_left:end_right,phase_unwrap_v(start_left:end_right),1);
    phase_unwrap_fit_v = phase_unwrap_v*0;
    phase_unwrap_fit_v(start_left:end_right) = ...
        (p(2)+p(1)*(start_left:end_right));
    %** Check if unwrap(angle(exp(1i*phase_unwrap_fit_v))) is the same as above.
    %** IT IS SAME.
    %phase_unwrap_fit_v(start_left:end_right) = ...
    %    unwrap(angle(exp(1i* (p(1)*(start_left:end_right)+p(2)) )));
    %figure,plot(phase_unwrap_fit_v)
    
    
    %* Expand fit region to 1:Nx
    x_v = 1:Nx;
    %x1 = start_left;
    %y1 = phase_unwrap_fit_v(x1);
    %x2 = end_right;
    %y2 = phase_unwrap_fit_v(x2);
    %y = (y2-y1)/(x2-x1)*(x_v-x1)+y1;
    y = polyval(p,x_v);
    phase_v = exp(1i*y);
    %** Test plot.
    %figure,plot(x_v,phase_unwrap_fit_v,'-r',x_v,y,':b')
    %figure,plot(x_v,angle(pi_subtr_avg_cut_v),'-r',...
    %  x_v,angle(exp(1i*y)),':b')
    %title(sprintf('pi\\_subtr\\_avg\\_cut\\_v(r), phase\\_v(b)'))
    
    
    %* Final phase correction matrix.
    %** Phase must be multiplied.
    %phase_mult_m = pi_subtr_m;
    phase_mult_m = repmat(phase_v,size(pi_subtr_m,1),1);
    
    
    %* Get out of this loop.
    flag_temp = 0;
end



%==========================================
% Fit to individual ky-line of phase.
% flag_fitPhase = 2.
while flag_temp==2
    
    phase_mult_m = pi_subtr_m * 0;
    
    for ind_ky = 1:length(ind_ind_ky_mult)
        
        % Take phase for each ky line.
        pi_subtr_each_v = pi_subtr_m(ind_ky,:);
        
        % Take only central phase values from pi_subtr_avg_v.
        center_length = round(length(pi_subtr_each_v)/4);
        take_center = floor(length(pi_subtr_each_v)/2)+1;
        start_left = take_center - floor(center_length/2);
        end_right = take_center + floor(center_length/2)-1;
        pi_subtr_each_cut_v = pi_subtr_each_v*0;
        pi_subtr_each_cut_v(start_left:end_right) = ...
            pi_subtr_each_v(start_left:end_right);
        
        % Fit the central part of phase to linear function.
        phase_unwrap_v = unwrap(angle(pi_subtr_each_cut_v)); % unwrap before fitting
        p = polyfit(start_left:end_right, ...
            phase_unwrap_v(start_left:end_right),1);
        %phase_unwrap_fit_v = phase_unwrap_v*0;
        %phase_unwrap_fit_v(start_left:end_right) = ...
        %    (p(2)+p(1)*(start_left:end_right));
        %* Check if unwrap(angle(exp(1i*phase_unwrap_fit_v))) is the same as above.
        %* THIS IS THE SAME.
        %phase_unwrap_fit_v(start_left:end_right) = ...
        %    unwrap(angle(exp(1i* (p(1)*(start_left:end_right)+p(2)) )));
        %figure,plot(phase_unwrap_fit_v)
        
        
        % Expand fit region to 1:Nx.
        x_v = 1:Nx;
        %x1 = start_left;
        %y1 = phase_unwrap_fit_v(x1);
        %x2 = end_right;
        %y2 = phase_unwrap_fit_v(x2);
        %y = (y2-y1)/(x2-x1)*(x_v-x1)+y1;
        y = polyval(p,x_v);
        phase_v = exp(1i*y);
        %** Test plot.
        %figure,plot(x_v,phase_unwrap_fit_v,'-r',x_v,y,':b')
        %figure
        %plot(x_v,angle(pi_subtr_each_cut_v),'-r',...
        % x_v,angle(exp(1i*y)),':b')
        %title(sprintf('pi\\_subtr\\_each\\_cut\\_v(r), phase\\_v(b), ky[%d]',ind_ky))
        
        
        % Reserve fit phase for each ky line.
        phase_mult_m(ind_ky,:) = phase_v;
    end
    
    
    % Final phase correction matrix.
    %* Phase must be multiplied.
    %phase_mult_m = pi_subtr_m;
    phase_mult_m = repmat(phase_v,size(pi_subtr_m,1),1);
    
    
    % Get out of this loop.
    flag_temp = 0;
end % while flag_temp



%==========================================
% Fit to individual ky-line of phase.
% Consider zeroth order phase of odd ky
% lines and its difference to even ky lines.
% flag_fitPhase = 3.
% DON'T USE THIS.
while flag_temp==3
    
    phase_mult_m = pi_subtr_m * 0;
    
    for ind_ky = 1:length(ind_ind_ky_mult)
        
        %* Take phase for each ky line.
        pi_subtr_each_v = pi_subtr_m(ind_ky,:);
        
        %* Take only central phase values from pi_subtr_avg_v.
        center_length = round(length(pi_subtr_each_v)/4);
        take_center = floor(length(pi_subtr_each_v)/2)+1;
        start_left = take_center - floor(center_length/2);
        end_right = take_center + floor(center_length/2)-1;
        pi_subtr_each_cut_v = pi_subtr_each_v*0;
        pi_subtr_each_cut_v(start_left:end_right) = ...
            pi_subtr_each_v(start_left:end_right);
        
        %* Take only central phase values from pi_odd_m.
        pi_odd_each_v = pi_odd_m(ind_ky,:);
        pi_odd_each_cut_v = pi_odd_each_v*0;
        pi_odd_each_v(start_left:end_right) = ...
            pi_odd_each_v(start_left:end_right);
        
        
        %* Fit the central part of phase to linear function.
        phase_unwrap_v = unwrap(angle(pi_subtr_each_cut_v)); % unwrap before fitting
        p = polyfit(start_left:end_right, ...
            phase_unwrap_v(start_left:end_right),1);
        phase_unwrap_fit_v = phase_unwrap_v*0;
        phase_unwrap_fit_v(start_left:end_right) = ...
            (p(2)+p(1)*(start_left:end_right));
        %** Check if unwrap(angle(exp(1i*phase_unwrap_fit_v))) is the same as above.
        %** IT IS SAME.
        %phase_unwrap_fit_v(start_left:end_right) = ...
        %    unwrap(angle(exp(1i* (p(1)*(start_left:end_right)+p(2)) )));
        %figure,plot(phase_unwrap_fit_v)
        
        %* Fit the central part of phase to linear function.
        phase_odd_unwrap_v = unwrap(angle(pi_odd_each_v)); % unwrap before fitting
        p_odd = polyfit(start_left:end_right, ...
            phase_odd_unwrap_v(start_left:end_right),1);
        
        %* Calculate zeroth order fit phase difference.
        p0_delta = mod(p_odd(2) - p(2),pi);
        
        
        %* Expand fit region to 1:Nx
        x_v = 1:Nx;
        %x1 = start_left;
        %y1 = phase_unwrap_fit_v(x1);
        %x2 = end_right;
        %y2 = phase_unwrap_fit_v(x2);
        %y = (y2-y1)/(x2-x1)*(x_v-x1)+y1;
        y = polyval(p,x_v);
        %phase_v = exp(1i*y);
        phase_v = exp(1i*y) * exp(1i*p0_delta);
        %** Test plot.
        %figure,plot(x_v,phase_unwrap_fit_v,'-r',x_v,y,':b')
        %figure
        %plot(x_v,angle(pi_subtr_each_cut_v),'-r',...
        % x_v,angle(exp(1i*y)),':b')
        %title(sprintf('pi\\_subtr\\_each\\_cut\\_v(r), phase\\_v(b), ky[%d]',ind_ky))
        
        
        %* Reserve fit phase for each ky line.
        phase_mult_m(ind_ky,:) = phase_v;
    end
    
    
    %* Final phase correction matrix.
    %** Phase must be multiplied.
    %phase_mult_m = pi_subtr_m;
    phase_mult_m = repmat(phase_v,size(pi_subtr_m,1),1);
    
    
    %* Get out of this loop.
    flag_temp = 0;
end % while flag_temp



%==========================================
% No fit. Use the whole pi_subtr_m.
% flag_fitPhase = 4.
while flag_temp==4
    
    %* Average magnitude data over ky direction.
    pi_odd_v = sum(abs(pi_odd_m),1);
    
    
    %* Take only portions of pi_odd_v to define object region.
    ind = find(abs(pi_odd_v) >= 0.1 * mean(abs(pi_odd_v)));
    object_v = pi_odd_v*0;
    object_v(ind) = pi_odd_v(ind);
    %** Test plot.
    %figure,plot(1:Nx,abs(pi_odd_v),'-r',1:Nx,abs(object_v),':b')
    %title(sprintf('pi\\_odd\\_v(r), object\\_v(b)'))
    %hold on
    %plot(ind(1),object_v(ind(1)),'or',ind(end),object_v(ind(end)),'ob')
    %hold off
    
    
    %* Expand above ind as much as 5% of total Nx.
    %expand_vox = ceil(0.05*Nx/2);
    %ind_expand = ind(1)-expand_vox:ind(end)+expand_vox;
    
    
    %* Don't use expanded index.
    ind_expand = ind;
    
    
    %--- Method 1.
    %* Take ind_expand portion of phase only from pi_subtr_m.
    phase_mult_m = pi_subtr_m * 0;
    phase_mult_m(:,ind_expand) = pi_subtr_m(:,ind_expand);
    %* Test plot.
    %showimage(jet,angle(phase_mult_m),angle(pi_subtr_m))
    
    
    %--- Method 2.
    %* Take average of phase_mult_m.
    %phase_mult_m = pi_subtr_m * 0;
    %phase_mult_m(:,ind_expand) = pi_subtr_m(:,ind_expand);
    %phase_v = mean(phase_mult_m,1);
    %phase_mult_m = repmat(phase_v,[size(phase_mult_m,1),1]);
    %* Test plot.
    %showimage(jet,angle(phase_mult_m),angle(pi_subtr_m))
    
    
    %* Get out of this loop.
    flag_temp = 0;
end



%% Phase correction

% Take the phase of 'phase_mult_m'.
% if strcmpi(DWIparams.shot_mode,'single-shot') && ...
%         ind_dw==0 && ind_avg==2
%     phase_mult_m = conj(phase_mult_m);
% end
    
    
% Generate ky index for +ve and -ve gradient for each shot in STD.
ky_shot_odd_ind_m = [];
ky_shot_even_ind_m = [];
for ind_shot = 1:DWIparams.nSHOT
    if ind_ec==1
        eval(sprintf('ky_shot%d_v = ind_shot:DWIparams.nSHOT:ny;',ind_shot))
    elseif ind_ec==2
        eval(sprintf('ky_shot%d_v = ind_shot:DWIparams.nSHOT:ny;',ind_shot))
    end
    eval(sprintf('ky_shot_odd_ind_m = [ky_shot_odd_ind_m; ky_shot%d_v(1:2:end)];', ...
        ind_shot))  % odd ky lines in each shot
    eval(sprintf('ky_shot_even_ind_m = [ky_shot_even_ind_m; ky_shot%d_v(2:2:end)];', ...
        ind_shot))  % even ky lines in each shot
end

% Do phase correction for each shot.
for ind_shot = 1:DWIparams.nSHOT
    k_corr_odd_temp_m = k_m(ky_shot_odd_ind_m(ind_shot,:),:);
    k_corr_even_temp_m = k_m(ky_shot_even_ind_m(ind_shot,:),:);
    
    h_odd_m = ift1(k_corr_odd_temp_m,2);
    h_even_m = ift1(k_corr_even_temp_m,2);
    
    
    % Use temp hybrid data in correction.
    h_odd_temp_m = h_odd_m;
    h_even_temp_m = h_even_m;
    
    
    % Match gradient polarity between averages.
    %* For SSH DTI, only b=0 has two averages. In this case, STD_start_sign
    %* is the same as PHC_start_sign for ind_avg==1, but different for 
    %* ind_avg==2.
    %* For MSH DTI, ind_avg increases up to the NSA. In this case, if 
    %* STD_start_sign is the same as PHC_start_sign, then it is the same
    %* for every ind_avg. Then compare STD_start_sign and PHC_start_sign rather
    %* considering ind_avg==1 or 2.
    if DWIparams.STD_start_sign==DWIparams.PHC_start_sign
        h_odd_corr_m = h_odd_temp_m;
        h_even_corr_m = h_even_temp_m .* phase_mult_m;
    else
        h_odd_corr_m = h_odd_temp_m;
        h_even_corr_m = h_even_temp_m .* conj(phase_mult_m);
    end
    
    
    % FT phase corrected data.
    k_odd_m = ft1(h_odd_corr_m,2);
    k_even_m = ft1(h_even_corr_m,2);
    
    
    % Take phase corrected data.
    k_corr_m(ky_shot_odd_ind_m(ind_shot,:),:) = k_odd_m;
    k_corr_m(ky_shot_even_ind_m(ind_shot,:),:) = k_even_m;
end

% showimage(jet,abs(k_corr_m),abs(ift2(k_corr_m)), ...
%    sprintf('epiCorrMethod=5,phx'))
% showimage(jet,abs(ift2(k_m)),angle(k_m))
% showimage(jet,abs(ift2(k_corr_m)),angle(ift2(k_corr_m)))%, ...
% sprintf('epiCorrMethod[%d],COIL[%d],fitPhase[%d]', ...
% DATAFLAGparams.epiCorrMethod,ind_coil,flag_fitPhase))



%% END






