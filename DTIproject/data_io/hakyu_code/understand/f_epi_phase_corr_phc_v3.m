
function k_corr_m = f_epi_phase_corr_phc_v3(k_m,phc_data,DATAFLAGparams, ...
    DWIparams,ind_sl,ind_chunk,ind_ec,ind_coil,ind_avg,ind_dw)
%[f_epi_phase_corr_phc_v3] does EPI echo phase correction using PHC data.
%
% USAGE
%   k_corr_m = f_epi_phase_corr_phc_v3(k_m,phc_data,DATAFLAGparams, ...
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
% 2011.10.19.
%   Commented 'if strcmpi(DWIparams.shot_mode,'single-shot')...' part for
%   rsEPI case. This may need to be uncommented in the future.
% 2012.02.01.
%   Take care when scan mode==3D.
% 2012.02.03.
%   Use different routine for finding pk_odd_m and pk_even_m.
% 2012.03.22.
%   Take care of single-shot with partial Fourier acqusition (Halfscan)
%   and/or NSA=1 for ind_avg==2 in b=0 data.
% 2012.05.16.
%   Add total variation routine to select the best EPI ghost corrected
%   image using the methods 1,2, and 4 (not 3). This is implemented in this
%   [f_epi_phase_corr_phc_v2.m].
% 2012.06.16.
%   Take care of scan_mode = 3D by taking ind_chunk as one of input
%   arguments.
%
% Ha-Kyu



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
if strcmpi(DWIparams.shot_mode,'multishot')
    ny = DWIparams.nSHOT * DWIparams.EPI_FACTOR(ind_ec);
else
    if strcmpi(DWIparams.halfscan,'yes')
        ny = size(k_m,1);
    end
end

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
    if strcmpi(DWIparams.scan_mode,'3D')
        p1_m = squeeze(phc_data(:,idx_ky,idx_kz,ind_chunk, ...
            idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,aver,1,:)).';
        p2_m = squeeze(phc_data(:,idx_ky,idx_kz,ind_chunk, ...
            idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,aver,2,:)).';
    else        
        p1_m = squeeze(phc_data(:,idx_ky,idx_kz,ind_sl, ...
            idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,aver,1,:)).';
        p2_m = squeeze(phc_data(:,idx_ky,idx_kz,ind_sl, ...
            idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,aver,2,:)).';
    end
elseif aver==3
    if strcmpi(DWIparams.scan_mode,'3D')
        p1_m = squeeze(mean(phc_data(:,idx_ky,idx_kz,ind_chunk, ...
            idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,:,1,:),10)).';
        p2_m = squeeze(mean(phc_data(:,idx_ky,idx_kz,ind_chunk, ...
            idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,:,2,:),10)).';
    else
        p1_m = squeeze(mean(phc_data(:,idx_ky,idx_kz,ind_sl, ...
            idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,:,1,:),10)).';
        p2_m = squeeze(mean(phc_data(:,idx_ky,idx_kz,ind_sl, ...
            idx_dyn,idx_card,ind_ec,idx_mix,ind_coil,:,2,:),10)).';
    end
    
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
if strcmpi(DWIparams.shot_mode,'multishot') % asumes two shot?????? fwb
    ind_odd_ky = 1:2:DWIparams.EPI_FACTOR(ind_ec);
    ind_even_ky = 2:2:DWIparams.EPI_FACTOR(ind_ec);
else % single-shot
    if strcmpi(DWIparams.halfscan,'yes')
        ind_odd_ky = 1:2:size(phc_data,12);
        ind_even_ky = 2:2:size(phc_data,12);
    end
end


% Read PHC data in +ve and -ve gradient.
%* This will make even ky lines are always smaller or equal to odd ky lines.
% if DWIparams.PHC_start_sign==1 && ~strcmpi(DWIparams.scan_mode,'3D')
%     pk_odd_m = p1_m(ind_odd_ky,:);
%     pk_even_m = p2_m(ind_even_ky,:);
% else
%     pk_odd_m = p2_m(ind_odd_ky,:);
%     pk_even_m = p1_m(ind_even_ky,:);
% end

%* Use different routine.
v = p1_m(1,:);
if all(v==0) % if top row is all-zero (b/c not collected?) - fwb
    pk_odd_m = p2_m(ind_odd_ky,:);
    pk_even_m = p1_m(ind_even_ky,:);
else
    pk_odd_m = p1_m(ind_odd_ky,:);
    pk_even_m = p2_m(ind_even_ky,:);
end

% p1/p2_m must differ by the direction of the gradient - fwb
% column direction is read-out (normal convention) - fwb

% Get hybrid-space data.
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
% if strcmpi(DWIparams.shot_mode,'single-shot') && ...
%         ind_dw==0 && ind_avg==2
%     pi_subtr_m = pi_even_m .* conj(pi_odd_m(ind_ind_ky_mult,:));
% end


% Normalize and make NaN and Inf zero.
pi_subtr_m = pi_subtr_m./abs(pi_subtr_m);
pi_subtr_m(isnan(pi_subtr_m)) = 0;
pi_subtr_m(isinf(pi_subtr_m)) = 0;



%% ==========================================
% Fit to ky-directionally averaged phase.
% flag_fitPhase = 1.
%while flag_temp==1
    
    %* Average phase over ky direction.
    pi_subtr_avg_v = sum(pi_subtr_m,1)/size(pi_subtr_m,1);
    %phase_v = pi_subtr_avg_v; % use phase_v after fitting to linear
    
    %* Average phase.
    %pi_subtr_avg_v = exp(1i*sum(angle(pi_subtr_m),1)/size(pi_subtr_m,1));
    
    
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
    %    1:Nx,angle(pi_subtr_avg_cut_v),':b')
    %title(sprintf('pi\\_subtr\\_avg\\_v(r), pi\\_subtr\\_avg\\_cut\\_v(b)'))
    
    
    %* Fit the central part of phase to linear function.
    phase_unwrap_v = unwrap(angle(pi_subtr_avg_cut_v)); % unwrap before fitting
    
    
    %** TEMP START ----------------------------------
%     v = del2(phase_unwrap_v);
%     ind = find(abs(v) > 0 & abs(v) < 0.05);
%     vv = phase_unwrap_v*nan;
%     vv(ind) = phase_unwrap_v(ind);
%     %figure,plot(phase_unwrap_v,'r'),hold on,plot(vv,'.b'),hold off
%     ind_v = find(~isnan(vv));
%     start_left = ind_v(1);
%     end_right = ind_v(end);
    %** TEMP END ----------------------------------
    
    p = polyfit(start_left:end_right,phase_unwrap_v(start_left:end_right),1);
    phase_unwrap_fit_v = phase_unwrap_v*0;
    phase_unwrap_fit_v(start_left:end_right) = (p(2)+p(1)*(start_left:end_right));
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
    % x_v,angle(exp(1i*y)),':b')
    %title(sprintf('pi\\_subtr\\_avg\\_cut\\_v(r), phase\\_v(b)'))
    
    
    %* Final phase correction matrix.
    %** Phase must be multiplied.
    %phase_mult_m = pi_subtr_m;
    phase_mult_m = repmat(phase_v,size(pi_subtr_m,1),1);
    
    
    %* Get out of this loop.
    flag_temp = 0;
%end
phase_mult1_m = phase_mult_m;



%% ==========================================
% Fit to individual ky-line of phase.
% flag_fitPhase = 2.
%while flag_temp==2
    
    phase_mult_m = pi_subtr_m * 0;
    
    for ind_ky = 1:length(ind_ind_ky_mult)
        
        % Take phase for each ky line.
        pi_subtr_each_v = pi_subtr_m(ind_ky,:);
        
        % Take only central phase values from pi_subtr_avg_v.
        %center_length = round(length(pi_subtr_each_v)/4);
        center_length = round(length(pi_subtr_each_v)/8);
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
    phase_mult_m = repmat(phase_v,size(pi_subtr_m,1),1); %looks like he throws away all the work from the loop here - fwb
    
   
    
    
    % Get out of this loop.
    flag_temp = 0;
%end % while flag_temp
phase_mult2_m = phase_mult_m;



%% ==========================================
% Fit to individual ky-line of phase.
% Consider zeroth order phase of odd ky
% lines and its difference to even ky lines.
% flag_fitPhase = 3.
% DON'T USE THIS.
%while flag_temp==3
    
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
%end % while flag_temp
phase_mult3_m = phase_mult_m;



%% ==========================================
% No fit. Use the whole pi_subtr_m.
% flag_fitPhase = 4.
%while flag_temp==4
    
    %* Average magnitude data over ky direction.
    pi_odd_v = sum(abs(pi_odd_m),1)/size(pi_odd_m,2); %pi_odd_m is the 
    
    
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
%end
phase_mult4_m = phase_mult_m;



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
    elseif ind_ec==2 % why??? - fwb
        eval(sprintf('ky_shot%d_v = ind_shot:DWIparams.nSHOT:ny;',ind_shot))
    end
    % append the matrices created above to the odd/even ind_m's - fwb
    eval(sprintf('ky_shot_odd_ind_m = [ky_shot_odd_ind_m; ky_shot%d_v(1:2:end)];', ...
        ind_shot))  % odd ky lines in each shot
    eval(sprintf('ky_shot_even_ind_m = [ky_shot_even_ind_m; ky_shot%d_v(2:2:end)];', ...
        ind_shot))  % even ky lines in each shot
end
% the final result is a pair of indexing matricies that tell the user the
% indicies of the odd/even lines (across cols) per shot (down rows). I think
% that the rows are all identical....??? - fwb

% Do phase correction for each shot: Use total variation in this version.
tot_var_v = zeros(1,4);
%lap_v = zeros(1,4); % Laplacian
%l2norm_v = zeros(1,4); % L2 norm
if DATAFLAGparams.epiCorrSelTotalVar==0
    method_v = flag_fitPhase;
else
    method_v = 0:4; % 0 is for no correction
end

% pick a method (they are all performed in the code above) - fwb
for ind_method = method_v   
    switch ind_method
        case 1
            phase_mult_m = phase_mult1_m;
        case 2
            phase_mult_m = phase_mult2_m;
        case 3
            phase_mult_m = phase_mult3_m;
        case 4
            phase_mult_m = phase_mult4_m;
        case 0
            phase_mult_m = ones(size(phase_mult1_m));
        otherwise
            error('Unknown ind_method')
    end
    
    % In case no correction is required.
    if flag_fitPhase==0
        phase_mult_m = ones(size(phase_mult_m));
    end
    
    % Phase correction.
    for ind_shot = 1:DWIparams.nSHOT
        k_corr_odd_temp_m = k_m(ky_shot_odd_ind_m(ind_shot,:),:);
        k_corr_even_temp_m = k_m(ky_shot_even_ind_m(ind_shot,:),:);
        
        % into image space - fwb
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
        
        std_start_sign = DWIparams.STD_start_sign;
        phc_start_sign = DWIparams.PHC_start_sign;
        
        %* Consider single-shot, NSA==1 case.
        if strcmpi(DWIparams.shot_mode,'single-shot') && ...
                (DWIparams.nNSA==1) && ind_avg==2 && ind_dw==0
            std_start_sign = -std_start_sign;
        end
        
        if std_start_sign==phc_start_sign
            h_odd_corr_m = h_odd_temp_m;
            h_even_corr_m = h_even_temp_m .* phase_mult_m;
        else
            h_odd_corr_m = h_odd_temp_m;
            h_even_corr_m = h_even_temp_m .* conj(phase_mult_m);
        end
        
        
        % FT phase corrected data. - back to kspace - fwb
        k_odd_m = ft1(h_odd_corr_m,2);
        k_even_m = ft1(h_even_corr_m,2);
        
        
        % Take phase corrected data. - make one image - fwb
        k_corr_m(ky_shot_odd_ind_m(ind_shot,:),:) = k_odd_m;
        k_corr_m(ky_shot_even_ind_m(ind_shot,:),:) = k_even_m;
    end
    
    eval(sprintf('k_corr%d_m = k_corr_m;',ind_method))
    
    % Calculate total variation on central region of aliased image.
    if DATAFLAGparams.epiCorrSelTotalVar==1
        i_corr_m = abs(ift2(k_corr_m));
        a = floor(size(k_corr_m,1)/3); %1/3 * nrows
        i_corr_m = i_corr_m((a+1):2*a,:);
        tot_var_v(ind_method) = f_total_var(i_corr_m);
        %lap_v(ind_method) = f_lap(i_corr_m);
        %l2norm_v(ind_method) = norm(i_corr_m(:));
    end
end

% Find minimum total variation or use default.
if DATAFLAGparams.epiCorrSelTotalVar==1
    ind = find(min(tot_var_v)==tot_var_v);
    if length(ind)>2
        ind = ind(1);
    end
else
    ind = flag_fitPhase; % use DATAFLAGparams.epiCorrFitMethod
end
eval(sprintf('k_corr_m = k_corr%d_m;',ind))

%fprintf('Min total variation is on method %d\n',ind)
%disp(tot_var_v)
%disp(lap_v)
%disp(lap_v.*tot_var_v)
%disp(l2norm_v+tot_var_v)

% showimage(jet,abs(k_corr_m),abs(ift2(k_corr_m)), ...
%    sprintf('epiCorrMethod=5,phx'))
% showimage(jet,abs(ift2(k_m)),angle(k_m))
% showimage(jet,abs(ift2(k_corr_m)),angle(ift2(k_corr_m)))%, ...
% sprintf('epiCorrMethod[%d],COIL[%d],fitPhase[%d]', ...
% DATAFLAGparams.epiCorrMethod,ind_coil,flag_fitPhase))



%% END






