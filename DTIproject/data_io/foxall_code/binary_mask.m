%
%  BINARY MASK
%
%  Provides a binary mask for image magnitude signal for either 2D or 3D
%  images.
% 
%  The method is based on Parseval's Theorem that states that sums of the 
%  squares of all the image pixels is equal to the sums of the squares of
%  all its kspace pixels. 
% 
%  The inverse FFT used to create the image kspace groups all the signal
%  into a tight peak located at the center of the kspace array. Summing  
%  kspace power over the whole array and just over its central region gives
%  estimates of the total power and the signal power. Noise power is 
%  estimated as the difference between the two values. Parseval's Theorem
%  and the additvity property of the Fourier transform allow these power
%  estimates for the noise field function and the signal field function to 
%  be applied in the image domain.
%
%  The pixels of the power image (the image squared) are sorted according
%  to size and an index lookup table for the sort is created to allow
%  mapping between the sorted and non sorted power lists. The cumulative 
%  power curve is computed for the sorted array.
%
%  Knowledge of the cumulative image power in terms of pixel power size is
%  then used to remove the set of lowest power pixels whose sum is equal 
%  to the noise power estimate made from the image kspace. The index in the   
%  look up table is then used to set the corresponding mask pixels equal 
%  to zero. 
% 
%  FORMAT 
% 
%  mask = binary_mask(im,threshold);
%
%  INPUT
%
%  im         float[]      Image array in either complex or magnitude form.
%
%  threshold  float        Remove a threshold factor (nominal = 1.0) times the 
%                          noise power from the cumulative distribution. 
%
%  OUTPUT
%
%  mask       int[]        Mask array of same dimension as the image having
%                          a value 0 for noise pixels or 1 for signal pixels.
%
%  $Id: binary_mask.m,v 1.1 2007/01/08 20:07:07 dfoxall Exp $
%

function mask = binary_mask(im,threshold);

%----------------------------------------------------------------
%  GENERAL OPERATIONS
%  * Compute the size and dimesnion of the input image array
%  * Compute the power image
%  * Compute the kspace for the magntiude image. 
%  * Compute the kspace power renormalizing the FFT with the
%    total number of image pixels.
%----------------------------------------------------------------

if (nargin < 2)
   threshold = 1.0;
end

im_SIZE   = size(im);
im_DIM    = max(size(im_SIZE));
im_PROD   = prod(im_SIZE);    
abs_im    = abs(im);
pwr_im    = abs_im.*abs_im;
ks        = fftshift(ifftn(fftshift(abs_im)));
abs_ks    = abs(ks);
pwr_ks    = im_PROD*abs_ks.*abs_ks;
mask      = ones(im_SIZE);


%----------------------------------------------------------------
% 2D METHOD
% * Sum the total kspace power.
% * Sum the kspace power over a central region with half the
%   image dimensions.
% * Estimate the noise power as the difference between the
%   two summations.
% * Make a 1D list of individual pixel power from the power image.
% * Sort the pixel power list by size and make a pixel index look
%   up table.
% * Compute the cumlative image power versus pixel power curve.
% * Remove those pixels from the mask for which the cumulative
%   image power is below the estimated noise power.  
%----------------------------------------------------------------
if ( im_DIM == 2 )
     beg_1  = round(im_SIZE(1)/4);
     beg_2  = round(im_SIZE(2)/4);
     fin_1  = round(3*im_SIZE(1)/4);
     fin_2  = round(3*im_SIZE(2)/4);
     pwr_t  = sum(sum(pwr_ks));
     pwr_s  = sum(sum(pwr_ks(beg_1:fin_1,beg_2:fin_2)));
     pwr_n  = (pwr_t - pwr_s);
   
     for I = 1:im_SIZE(1)
     for J = 1:im_SIZE(2)
         N            = I + im_SIZE(1)*(J-1);
         pwr_list(N)  = pwr_im(I,J);
     end
     end

     [pwr_sort,pwr_ind] = sort(pwr_list);
     cum_pwr            = cumsum(pwr_sort);
     pwr_lim            = pwr_n*ones(size(pwr_sort));

     for N = 1:im_PROD
         if ( cum_pwr(N) <= threshold*pwr_n )
            ind_N = pwr_ind(N);
            ind_i = 1 + rem(ind_N-1,im_SIZE(1));
            ind_j = 1 + round((ind_N - ind_i)/im_SIZE(1));
            mask(ind_i,ind_j) = 0;
         end 
     end
end

%----------------------------------------------------------------
% 3D METHOD
% Simple extension of the 2D method allowing for the extra
% extra dimension.
%----------------------------------------------------------------
if ( im_DIM == 3 )
     beg_1  = round(im_SIZE(1)/4);
     beg_2  = round(im_SIZE(2)/4);
     beg_3  = round(im_SIZE(3)/4);
     fin_1  = round(3*im_SIZE(1)/4);
     fin_2  = round(3*im_SIZE(2)/4);
     fin_3  = round(3*im_SIZE(3)/4);
     pwr_t  = sum(sum(sum(pwr_ks)));
     pwr_s  = sum(sum(sum(pwr_ks(beg_1:fin_1,beg_2:fin_2,beg_3:fin_3))));
     pwr_n  = (pwr_t - pwr_s);

     for I = 1:im_SIZE(1)
     for J = 1:im_SIZE(2)
     for K = 1:im_SIZE(3)
         N           = I + im_SIZE(1)*(J-1) + im_SIZE(1)*im_SIZE(2)*(K-1);
         pwr_list(N) = pwr_im(I,J,K);  
     end
     end
     end

     [pwr_sort,pwr_ind] = sort(pwr_list);
     cum_pwr            = cumsum(pwr_sort);
     pwr_lim            = pwr_n*ones(size(pwr_sort));

     for N = 1:im_PROD
         if ( cum_pwr(N) <= threshold*pwr_n )
            ind_N  = pwr_ind(N);
            ind_i  = 1 + rem(ind_N-1,im_SIZE(1));
            ind_j  = 1 + rem(round((ind_N - ind_i)/im_SIZE(1)),im_SIZE(2));
            ind_k  = 1 + round((ind_N - ind_i - (ind_j-1)*im_SIZE(1))/(im_SIZE(1)*im_SIZE(2)));
            mask(ind_i,ind_j,ind_k) = 0;
         end 
     end
end

%-----------
% END
%-----------
