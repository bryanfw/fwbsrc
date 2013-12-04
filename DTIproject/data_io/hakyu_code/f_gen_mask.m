
function mask_3d = f_gen_mask(t2_3d,const)
%[f_gen_mask] generates mask from T2 (b=0) image.
%
% Usage:
%   mask_3d = f_gen_mask(t2_3d,const)
%
% Note:
%   This is specially made for using 7 T non-DW image.
%
% See also f_gen_mask2  gen_mask
%
%
% Last modified
% 2011.01.11.
% 2011.01.26.
%   Add 2-D process which is faster than 3-D process.
% 2011.04.22.
%   Add const as an input parameter.
% 2011.12.18.
%   Take care for a < 0 error.
%
%
% Ha-Kyu


%% Main
% tic
% t2_3d = t2_3d/max(t2_3d(:));
% const = 0.1; % 0.1 will be the same as long as using DTI_NAV recon routine is used at 7 T
% mask = t2_3d > (graythresh(t2_3d)*const);
mask1 = abs(t2_3d)>0;
mask_3d = zeros(size(t2_3d),'single');
for ind_slice = 1:size(mask_3d,3)
    % 3-D.
    %mask = t2_3d > (graythresh(t2_3d)*const);
    
    % 2-D.
    t2_m = t2_3d(:,:,ind_slice);
    t2_m = t2_m/max(t2_m(:));
    
    % Stretch intensity - don't use.
    %lowhigh = stretchlim(t2_m);
    %t2_m = imadjust(t2_m,[lowhigh(1)/2,lowhigh(2)/2],[]);
    
    mask_m = t2_m > (graythresh(t2_m)*const);
    
    % Search mask.
    flag = 1;
    a = const;    
    while flag
        % 3-D.
        %m = mask(:,:,ind_slice);
        
        % 2-D.
        m = mask_m;        
        
        % Search.
        m = imfill(m,'holes');
        m1 = mask1(:,:,ind_slice);
        m1 = imerode(m1,strel('disk',9));
        mm = m & m1;
        L = bwlabel(mm);
        s = regionprops(L,'Area');
        c = struct2cell(s);
        v = [c{:}];
        ind = find(max(v)==v);
        if ~isempty(ind), ind = ind(1); end % in case there are more than one region
        if isempty(ind) || (v(ind) < 0.3*numel(m)) % 40 percent of image voxels 
            a = a-0.001; %[a, v(ind)]
            
            % 3-D.
            %mask = t2_3d > (graythresh(t2_3d)*a);
            
            % 2-D.
            mask_m = t2_m > (graythresh(t2_m)*a);
            
            % In case a < 0, though it will never happen.
            if a < 0
                %error('f_gen_mask:main','a < 0')
                a = 0;
            end
        else
            flag = 0;
        end
    end    
    m0 = (L==ind);
    mask_3d(:,:,ind_slice) = m0;
end
% fprintf('finished in %f sec\n',toc)
clear  t2_3d  mask  mask1 m m1 mm L s c v ind m0




%% END




