function slice_slider(im, mask,limits)
% SLICE_SLIDER(3d_volume,mask,range) plots a 3d volume 1 slice at a time
%   on the current figure.
% A slider bar allows you to slice through the volume.
% Assumes that you want to look at the third dimension using
%   imagesc and the default colormap.
% Fix it/make it better if you want to.
%
%
% Frederick Bryan, 2013

% 1st argument check
if nargin<1
    error('Not enough arguments.');
end
im = squeeze(im);
if ndims(im)>3 || ndims(im)<2
    error('data must be 2D or 3D');
end
if ~isreal(im);
    im = abs(im); 
    warning('modulus taken of imaginary data')
end

% 2nd and third arg checks
if nargin < 2;
    mask = [];  
end
if isempty(mask)
    %nothing, pass
else
    mask = squeeze(mask);
    if ~isequal(size(im), size(mask))
        error('mask and im must be same size');
    end
end
if nargin<3
    limits =[];
end

% build figure stuff
data.im = im;
clear im;

data.mid = ceil(size(data.im,3)/2);
data.top = size(data.im,3);
data.step = [1, 1] / (data.top - 1);
if any(isinf(data.step));
    data.step = [1,1];
end

data.figure1=gcf;
clf(data.figure1,'reset');

% hold on;

% create axes
% keyboard;
data.axes1 = axes;

% less than 400 unique cases and no decimal numbers
if ~isempty(limits)
    % use limits
elseif length(unique(data.im(:))) < 400 && ~isa(data.im,'logical') && ...
        isequal(round(data.im),data.im)
    limits = [min(data.im(:)) max(data.im(:))];
    colormap jet;
else
    if ~isa(data.im,'logical');
        data.im = img_normalize(data.im);
    end
    limits = [0 1];
    colormap gray;
    if ~isempty(mask);
        data.im = applymask(data.im,mask,limits);
    end
end

if limits(2) <= limits(1); % fix limits (all 0 special case)
    limits = [0 1];
end
% keyboard;
data.limits = limits;

data.ui = uicontrol('Style','slider','Min',1,'Max',data.top,'Value',data.mid,...
    'SliderStep', data.step, 'Position',[10 10 150 25],...
    'Callback',{@changeslice});

% update gui object
guidata(data.figure1,data);

changeslice(data.ui)

function changeslice(src,~)
% keyboard;
data = guidata(src);

% update slice by reploting figure
newsl = get(src,'Value');
axes(data.axes1);
imagesc(squeeze(data.im(:,:,round(newsl),:)),data.limits);
axis image; axis off;

% reprint title
title(sprintf('%3.0i',round(newsl)));

function im4 = applymask(im3, seg3,intlim)

im4 = zeros([size(im3) 3]);

max_label_num = max(seg3(:))+1;
im = double(im3);
mi = min(im3(:));
mx = max(im3(:));

% fix the intensity limits
intlim(intlim < 0) = 0;
intlim(intlim > 1) = 1;

cmap = jet;
alpha = .25;

for sl = 1:size(im3,3);
    
    im = im3(:, :, sl);
    seg = uint8(seg3(:, :, sl));
    
    
    % convert the raw image to an rgb image
    xl = mi + intlim(1) * (mx - mi);
    xh = mi + intlim(2) * (mx - mi);
    im(im < xl) = xl;
    im(im > xh) = xh;
    im = (im - xl) / (xh - xl);
    rgb = repmat(im, [1 1 3]);
    
    % resample the colormap
    cs = size(cmap, 1);
    inds = round(linspace(1, cs, max_label_num));
    cmap = cmap(inds, :);
    
    % convert the segmented image to an rgb image
    clr = double(ind2rgb(seg, cmap));
    
    % create the overlay image
    overlayim = rgb;
    inds = find(seg > 0);
    for ch = 1:3
        rgbch = rgb(:, :, ch);
        overlaych = overlayim(:, :, ch);
        clrch = clr(:, :, ch);
        overlaych(inds) = alpha * clrch(inds) + (1 - alpha) * rgbch(inds);
        overlayim(:, :, ch) = overlaych;
    end
    
    im4(:,:,sl,:) = overlayim;
    
end

function [out] =img_normalize(in,prctiles,mask)
% out = IMG_NORMALIZE(in)
% utility function used to normalize intensity to between 0.5 and 99.5
% percentile and remap everyting [0 1]

out = in;

if nargin < 2 || isempty(prctiles)
    prctiles = [.5 99.5];
end
if nargin==3 % if mask given, mask image
    if isa(mask,'logical') 
        in = in(mask);
    elseif all(ismember(mask(:),[0 1])); % check if all elements are in {0,1}
        mask = logical(mask);
        in = in(mask);
    else
        error('mask elements must be [0 1] if not class:logical');
    end
end

% cast to double
in = double(in);

% saturate btwen .5 and 99.5 percentile
lowval = prctile(in(:), prctiles(1));
highval = prctile(in(:), prctiles(2));

out(out<lowval) = lowval;
out(out>highval) = highval;

% scale to one;
out = (out-lowval)/(highval-lowval);

