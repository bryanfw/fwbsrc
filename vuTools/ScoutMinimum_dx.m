function scout_dx = ScoutMinimum_dx(kSpace,kStart,kStop,imgkStart,imgkStop,minMaxSearch)
% Course scout of image metric space to find minimum to initialize
% optimizer.
% Reference: "Improved optimization strategies for autofocusing motion
% compensation in MRI via the analysis of image metric maps." Wei Lin, Hee
% Kwon Song.
% By : Kevin Wilson
% Date : June 15, 2007
% Parameters:
%   kSpace : the k-space image (shifted)
%   kStart : the beginning line of k-space to be shifted
%   kStop : the end line of k-space to be shifted
%   imgkStart : the beginning line of k-space used to create the image
%   imgkStop : the end line of k-space used to create the image
%   minMaxSearch : search +/- this amount (fraction of image)

% delta x scout

% Calculate Resolution (one pixel)
dims = size(kSpace);
dt = (100/dims(2))/100;
% Store current min shift
currentMinT = [];
currentMinE = [];
% Calculate oscillation period Dy
Dy = 1/mean([kStart:kStop]);
% Loop over one oscillation period at 1/4 step size
for s = -Dy/2:Dy/4:Dy/2
    % Shift kSpace (in y)
    shift = exp(-i*2*pi/dims(1) .* [0:dims(1)-1] .* s * dims(1))';
    shift = repmat(shift,[1 dims(1)]);
    shift = shift(kStart:kStop,:);
    shifted_kSpace = kSpace;
    shifted_kSpace(kStart:kStop,:) = kSpace(kStart:kStop,:).*shift;

    % Loop over search space (in x)
    for t = -minMaxSearch:dt:minMaxSearch
        E = EntropyMetric(t,shifted_kSpace,kStart,kStop,imgkStart,imgkStop,1);
        if (isempty(currentMinE) || (E < currentMinE))
            currentMinE = E;
            currentMinT = t;
        end
    end
end

scout_dx = currentMinT;