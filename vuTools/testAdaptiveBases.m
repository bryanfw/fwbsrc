

gridResolutions = {[3 5],[5 8],[10 20]};
rbfRadiusPercent = [0.75 0.75 0.75];
gradientThreshPercent = [1 1 1];

im = imread('cameraman.tif');
%im = imresize(im,[245 206]);
im = vuGenerateMetaImage(im);
def = imread('cameraman_deformed.tif');
%def = imresize(def,[245 206]);
def = vuGenerateMetaImage(def);
mask = def;
mask.Data(:) = 0;
mask.Data(50:200,50:200)=1;

tic
[reg,XX,YY] = vuAdaptiveBasesRegistration(im,def, ...
    'gridResolutions',gridResolutions);
toc