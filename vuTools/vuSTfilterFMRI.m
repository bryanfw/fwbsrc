function img4D = vuSTfilterFMRI(im, ind_r, ind_i, fc, fs, verb)

% This function accepts an image for an fMRI experiment and returns
% the Stockwell transform filtered data as described in Goodyear et al.,
% Magnetic Resonance in Medicine, 51:16-21, 2004.  If the data are complex,
% then the magnitude and phase are filtered separately and recombined.
% 
% img4D = vuSTfilterFMRI(im, ind_r, ind_i, fc, fs, verb, sopt, spath);
% 
% INPUTS:
% 
% (1) im    - fMRI data
% (2) ind_r - Index in REC file for magnitude (or real) data (DEFAULT = 1)
% (3) ind_i - Index in REC file for imaginary data (DEFAULT = -1)
% (4) fc    - Scalar to define filter cut-off in Stockwell Transform domain
%             (DEFAULT = 3; see Goodyear et al., 2004) Frequency components
%             with a magnitude greater than fc times the median magnitude
%             for that frequency are replaced by fs * (median magnitude).
%             (Decreasing fc lowers the threshold to identify artifacts.)
% (5) fs    - Scalar to further suppress frequency components identified as
%             artifacts (DEFAULT = 1) For example, fs=0.9 would replace
%             these frequency components with 90% of the median magnitude.
% (6) verb  - Is verbose output required? (DEFAULT = 0)
% 
% OUTPUT:
% 
% img4D     - Image after Stockwell transform filtering (x, y, z, time).
%             Data will be complex if both real and imaginary components
%             are specified in the input.
% 
% NOTES:
% 
% As an example, consider a PAR/REC image containing all four forms of data:
% magnitude, real, imaginary, and phase -- in that order.  The command to
% process only the magnitude data is:
% 
%    im = vuOpenImage('myFile.PAR');
%    img4D = vuSTfilterFMRI(im, 1);
% 
% The command to process the complex data is:
% 
%    im = vuOpenImage('myFile.PAR');
%    img4D = vuSTfilterFMRI('myfile.PAR', 2, 3);
% 
% In this case, the magnitude and phase data are not required because they
% are calculated from the real and imaginary components.  If the filtered
% images are to be saved as a PAR/REC structure, then only indices 2 and 3
% are changed (i.e., indices 1 and 4 still contain the unfiltered images).
% 
% Finally, the command to verbosely process only the magnitude data using
% the default filter is:
%
%    im = vuOpenImage('myFile.PAR');
%    img4D = vuSTfilterFMRI('myfile.PAR', 1, -1, 3, 1, 1);
% 
% Copyright (c) 2008 - Vanderbilt University Institute of Imaging Science

% HISTORY
% 2008-08-15 - 1.02 - added 'fs' parameter (RLB)
% 2008-08-13 - 1.01 - able to output PAR/REC structure (RLB)
% 2008-08-08 - ver. 1.00 - implemented by Robert L. Barry

if nargin == 0
    disp_help;
    return;
end

if nargin < 6
    verb = 0;
    if nargin < 5
        fs = 1;
        if nargin < 4
            fc = 3;
            if nargin < 3   % assume data are not complex if
                ind_i = -1; % imaginary part is not specified
                if nargin < 2
                    ind_r = 1;
                end
            end
        end
    end
end
if ind_i < 1
    magnitude_only = 1;
else
    magnitude_only = 0; % imaginary part specified, so data are complex!
end

% It must be emphasized that this filter is applied to individual voxels
% in an fMRI data set (typically 80,000 for low-resolution and 200,000+ for
% high-resolution), so it is of the UTMOST IMPORTANCE that this code is as
% efficient as possible.

% Take only what we need and get rid of the rest
if magnitude_only
    img4D = squeeze(im.Data(:,:,:,ind_r,:));
else
    img4D = squeeze(im.Data(:,:,:,ind_r,:) + 1i * im.Data(:,:,:,ind_i,:));
end

spc = im.Spc;
origin = im.Origin;

clear im;

% Stockwell Transform parameters -- adapted from st.m --------------------
n = size(img4D,4);
maxfreq = fix(n/2);
vector = zeros(2,n);
vector(1,:) = 0:1:n-1;
vector(2,:) = -n:1:-1;
t2 = vector .* vector;
gw = zeros(maxfreq,n);
for v = 1:1:maxfreq % define Gaussian windows
    gw(v,:) = sum(exp(t2 * (-2*pi*pi/v/v)),1);
end
% ------------------------------------------------------------------------

arti = zeros(1,n);
PT_STF = zeros(n,1);
n1 = floor(n/2)+1; % determine where the Fourier spectrum reflects
n2 = n1 + mod(n,2);

for i1 = 1 : size(img4D,3) % process each slice separately
    if verb
        disp(sprintf('Processing slice %d of %d...', i1, size(img4D,3)));
    end
    for i2 = 1 : size(img4D,1)
        for i3 = 1 : size(img4D,2)

pt = double(squeeze(img4D(i2,i3,i1,:))); % single point (may be complex)

            for mp = 1 : 2-magnitude_only
                if mp == 1
                    PT = fft(abs(pt));   % process magnitude
                else
                    PT = fft(angle(pt)); % process phase
                end
        
% Stockwell Transform -- The following few lines (and the ST parameters
% defined above) were taken from st.m (written by R. G. Stockwell) and
% adapted for optimal efficiency.  (Refer to Stockwell et al., IEEE Trans.
% on Signal Processing, vol. 44, number 4, April 1996, pages 998-1001.)

pt_v = [PT, PT];
pt_ST = zeros(ceil(maxfreq+1), n);
if mp == 1, pt_ST(1,:) = mean(abs(pt));
else        pt_ST(1,:) = mean(angle(pt)); end
for loopy = 1:1:maxfreq
    pt_ST(loopy+1,:) = ifft(pt_v(loopy+1:loopy+n) .* gw(loopy,:));
end

% Identify artifacts for each frequency component.  As per Goodyear et al.
% (2004), an artifact begins when the magnitude exceeds three times the
% median and ends when the magnitude drops below the median.  These points
% are then replaced by the median magnitude of that frequency.

H_filtered = abs(pt_ST);     % we only need the magnitude (for now)
H_median = median(H_filtered,2); % calculate median for each frequency

for i4 = 1 : size(H_filtered,1) % examine one frequency component at a time
    med = H_median(i4);
    lim = fc * med;          % default filter cut-off is three times median
    Hi = H_filtered(i4,:);
    Ha = Hi > lim;           % label definite artifacts
    if any(Ha)               % can we skip this frequency?
        Ha = single(Ha);
        Ha1 = find(Ha == 1);
        if Ha1(end) == n, Ha1 = Ha1(1:end-1); end % ignore if last point
        ia = 0;              % reset index for number of artifacts
        for i5 = 1 : size(Ha1,2)
            if ~Ha(1,Ha1(1,i5)+1)
                ia = ia + 1;
                arti(1,ia) = Ha1(1,i5); % find end of artifacts
            end % if
        end % i5
        for i5 = 1 : ia
            i6 = arti(1,i5) + 1;
            while ((i6 <= n) && (Ha(i6) == 0) && (Hi(i6) > med))
                Ha(i6) = 1;  % label these points as artifacts too
                i6 = i6 + 1;
            end % while
        end % i5
        H_filtered(i4,:) = (Hi .* ~Ha) + ((fs*med) .* Ha); % apply filter
    end % skipfreq
end % i4

pt_ST_filtered = H_filtered .* exp(1i * angle(pt_ST));
PT_STF(1:n1,1) = sum(pt_ST_filtered, 2); % sum frequencies over time
PT_STF(n2:end,1) = conj(flipud(PT_STF(2:n1,1))); % Fourier transform
pt_stf = real(ifft(PT_STF)); % ST filtered profile

if mp == 1, mag_stf = pt_stf;
else        pha_stf = pt_stf; end

            end % mp

% Output filtered image
if magnitude_only
    img4D(i2,i3,i1,:) = reshape(single(mag_stf), [1 1 1 n]);
else
    img4D(i2,i3,i1,:) = reshape(single(mag_stf .* exp(1i * pha_stf)), ...
        [1 1 1 n]);
end
            
        end % i3
    end % i2
end % i1

img4D = vuGenerateMetaImage(img4D,spc,origin);
