close all
clear all

im = vuGenerateMetaImage(imread('cameraman.tif'));

ipp_im = MEXAllocateAlignedMemory(prod(im.Dims),0);
MEXCopyArrayToAlignedMemory(im.Data,ipp_im);

MEXDebugPrintAlignedMemory(ipp_im);

nBins = 32;
values = MEXAllocateAlignedMemory(nBins+1,0);
fixed_levels = MEXAllocateAlignedMemory(nBins+1,0);
moving_levels = MEXAllocateAlignedMemory(nBins+1,0);

MEXIPPCalculateImageBins(ipp_im,ipp_im,values,fixed_levels, ...
                         moving_levels,nBins);
fixed_entropy = MEXIPPCalculateFixedImageEntropy(ipp_im,values, ...
                                                 fixed_levels,nBins);

sample_in_pdf = MEXAllocateAlignedMemory(nBins*nBins,0);
sample_out_pdf = MEXAllocateAlignedMemory(nBins*nBins,0);

MEXIPPCalculateSampleJointPDF(ipp_im,ipp_im,values,fixed_levels, ...
                              moving_levels,sample_out_pdf,nBins);

MEXZeroAlignedMemory(sample_in_pdf)

moving_pdf = MEXAllocateAlignedMemory(nBins,1);
ln_joint_pdf = MEXAllocateAlignedMemory(nBins*nBins,0);
ln_moving_pdf = MEXAllocateAlignedMemory(nBins,1);
MEXIPPCalculateLimitedMemoryMutualInformation(fixed_entropy, sample_in_pdf, ...
                                           sample_out_pdf, moving_pdf, ...
                                           ln_joint_pdf, ln_moving_pdf, ...
                                           prod(im.Dims), 32);
MEXFreeAlignedMemory(ipp_im);
MEXFreeAlignedMemory(values);
MEXFreeAlignedMemory(fixed_levels);
MEXFreeAlignedMemory(moving_levels);
MEXFreeAlignedMemory(sample_in_pdf);
MEXFreeAlignedMemory(sample_out_pdf);
MEXFreeAlignedMemory(moving_pdf);
MEXFreeAlignedMemory(ln_joint_pdf);
MEXFreeAlignedMemory(ln_moving_pdf);