

gridResolutions = {[5]};
rbfRadiusPercent = [0.75];
gradientThreshPercent = [1];

im = imread('cameraman.tif');
im = imresize(im,0.25);
im(:,:,2) = im;
im(:,:,3:4) = im;
im(:,:,5:8) = im;
im(:,:,9:16) = im;
%im = imresize(im,[245 206]);
im = vuGenerateMetaImage(im);
def = imread('cameraman_deformed.tif');
def = imresize(def,0.25);
def(:,:,2) = def;
def(:,:,3:4) = def;
def(:,:,5:8) = def;
def(:,:,9:16) = def;
%def = imresize(def,[245 206]);
def = vuGenerateMetaImage(def);
mask = def;
mask.Data(:) = 0;
mask.Data(50:200,50:200,50:200)=1;

tic
[reg,XX,YY] = vuAdaptiveBasesRegistration(im,def, ...
    'gridResolutions',gridResolutions);
toc