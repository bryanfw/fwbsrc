function E = EntropyMetric(t,kSpace,kStart,kStop,imgkStart,imgkStop,direction)
% Metric for optimization
% Parameters:
%   t : time shift to be optimized (fraction of 1)
%   kSpace : the k-space image (shifted)
%   kStart : the beginning line of k-space to be shifted
%   kStop : the end line of k-space to be shifted
%   imgkStart : the beginning line of k-space used to create the image
%   imgkStop : the end line of k-space used to create the image
%   direction : the direction of the shift (x=1,y=2)

% Shift the correct lines in k-space (Shift Theorem)
dims = size(kSpace);
if (direction == 1)
    shift = exp(-i*2*pi/dims(2) .* [0:dims(2)-1] .* t * dims(2));
    shift = repmat(shift,[kStop-kStart+1 1]);
    kSpace(kStart:kStop,:) = kSpace(kStart:kStop,:).*shift;

    % Transform img and calculate entropy
    img = ifft2(kSpace(imgkStart:imgkStop,:));
    Bp = abs(img);
    % grad_img = vuGradientMagnitude(img);
    % Bp = abs(grad_img);
    Bp(Bp==0) = eps; % Correct bad behaving zeros
    Bmax = sqrt(sum(Bp(:).*Bp(:)));
    E = -1*sum((Bp(:)./Bmax).*log(Bp(:)./Bmax));
elseif (direction == 2)
    shift = exp(-i*2*pi/dims(1) .* [0:dims(1)-1] .* t * dims(1))';
    shift = repmat(shift,[1 dims(1)]);
    shift = shift(kStart:kStop,:);
    kSpace(kStart:kStop,:) = kSpace(kStart:kStop,:).*shift;

    % Transform img and calculate entropy
    img = ifft2(kSpace(:,imgkStart:imgkStop));
    Bp = abs(img);
    % grad_img = vuGradientMagnitude(img);
    % Bp = abs(grad_img);
    Bp(Bp==0) = eps; % Correct bad behaving zeros
    Bmax = sqrt(sum(Bp(:).*Bp(:)));
    E = -1*sum((Bp(:)./Bmax).*log(Bp(:)./Bmax));
else
    error('MATLAB:EntropyMetric','Unreconized direction');
end