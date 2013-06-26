function lab_vol = quick_label(invol)
% label_volume = QUICK_LABEL(intensity_vol)
% -or-
% label_volume_filename.nii.gz = QUICK_LABEL(input_file.nii/.nii.gz)
% Brings up the quick_label GUI, allowing you to rapidly manually label
% a volume with the help of logistic regression. 
% 
% If you pass a NIFTI file, it tries to do:
% axial slices
% anterior (nose-side of head) be up in the 2D slice
% anterior->superior (foot->head) be increasing slice number
% get dimensions correct
% 
% If you pass a regular matlab matrix, it assumes that x,y,z are left-right,up-down
%     and slice dimensions, respectively.
% 
% Instructions for use:
% GUI will pop up. 
% Navigate to the top of the anatomy (or one slice before) using the slider. 
% Outline the anatomy in question. 
% Continue for the rest of the anatomy. 
% When you've enclosed the anatomy completely (100% enclosed)
%     hit the Logisticize Me! button. 
%     We then try to classify the image based on the information you provided. 
% If you are happy with this classification, hit Finish! to return the logical 
%     variable lab_vol to the matlab workspace. 
% If not, modify the image by:
%     1. Clipping of stuff that shouldn't be included: start outside the included 
%        area and end outside.
%     2. Including stuff that should be included: Start inside the included area 
%        and end inside. 
%     Edge-cases: If a slice is totally wrong, "Redo-Slice" to draw a new polygon.
%     
% At any point, you can hit 3D render and see a 3D rendering of the current 
%     segmentation in a new figure window.
%     
% Notes:
% Homology is enforced.
% Nothing you exclude will ever be included. 


if isa(invol,'char')
    if isequal(invol(end-3:end),'.nii');
        vol = load_nii(invol);
    elseif isequal(invol(end-6:end),'.nii.gz')
        vol = load_nii_gz(invol);
    else
        fprintf('\nInput was string, but not .nii or .nii.gz. Exiting\n\n');
        return;
    end
    invol = vol.img; 
end

if isnumeric(invol)  % a matrix
    if ~(ndims(invol)>1 && ndims(invol)<4)
        fprintf('\nVolumes must be 2-3D. Exiting\n\n');
        return;
    end
else
    fprintf('\nVolume was non-numeric. Exiting\n\n');
    return;
end

keyboard;

mid = ceil(size(invol,3)/2);
top = size(invol,3);
step = [1, 1] / (top - 1);

h = figure; 

