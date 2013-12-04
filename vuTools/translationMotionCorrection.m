function corr_kSpace = translationMotionCorrection(kSpace,nLines,initialBaseWidth)
% Translation motion correction of kSpace data
% Based on "Automatic Compensation of Motion Artifacts in MRI"
% By: Kevin Wilson
% Date : June 15, 2007
% Parameters :
%   kSpace : Data to be corrected
%   nLines : number of lines to use as a segment during correction
%   initialBaseWidth : width of initial k-space to create base image

% Constants
dims = size(kSpace);
options = optimset('TolX',1e-8,'TolFun',1e-2,'DiffMaxChange',1e-2,'DiffMinChange',1e-8);
% middle of ky
midKy = floor(dims(1)/2);

% Variables
baseKSpaceRange =  [(midKy - floor(initialBaseWidth/2)) (midKy + floor(initialBaseWidth/2) - 1)];

% Flag to determine what side of base we're on
positiveSide = 1;

% Loop over image
while (baseKSpaceRange(1)>0) && (baseKSpaceRange(2)<dims(1))
    if (positiveSide)
        % Positive side
        segkStart = min(dims(1),baseKSpaceRange(2) + 1);
        segkStop = min(dims(1),segkStart + nLines - 1);
    else
        % Negative side
        segkStop = max(1,baseKSpaceRange(1) - 1);
        segkStart = max(1,segkStop - nLines + 1);
    end
    
    % Scout X
    if (positiveSide)
        scout_dx = ScoutMinimum_dx(kSpace,segkStart,segkStop,baseKSpaceRange(1),segkStop,0.1);
    else
        scout_dx = ScoutMinimum_dx(kSpace,segkStart,segkStop,segkStart,baseKSpaceRange(2),0.1);
    end
    
    % Loop though searches in dx and dy until no shift
    dx = 0.01;
    dy = 0.01;
    while(abs(dy)>0.0001&&abs(dx)>0.0001)
        dx = fminsearch(@(dx) EntropyMetric(dx,kSpace,segkStart,segkStop,baseKSpaceRange(1),baseKSpaceRange(2),1),scout_dx,options);
        
        % reset scout_dx, to start search at 0 next time
        scout_dx = 0;
        
        % Shift x axis
        shift = exp(-i*2*pi/256 .* [0:256-1] .* dx * 256);
        if (positiveSide)
            shift = repmat(shift,[dims(2)-segkStart+1 1]);
            kSpace(segkStart:dims(2),:) = kSpace(segkStart:dims(2),:).*shift;
        else
            shift = repmat(shift,[segkStop 1]);
            kSpace(1:segkStop,:) = kSpace(1:segkStop,:).*shift;
        end

        % Scout y at oscillations
        minE = [];
        Dy = 1/mean([segkStart:segkStop]);
        for a = 0:Dy:0.1
            % Search at i
            scout_dy = fminsearch(@(dy) EntropyMetric(dy,kSpace,segkStart,segkStop,baseKSpaceRange(1),baseKSpaceRange(2),2),a,options);
            % Find Entropy
            E = EntropyMetric(scout_dy,kSpace,segkStart,segkStop,baseKSpaceRange(1),baseKSpaceRange(2),2);
            if(isempty(minE)||E<minE)
                % store dy and E
                dy = scout_dy;
                minE = E;
            end
            % if i ~= 0 search negative side
            if (a~=0)
               % Search at -i
                scout_dy = fminsearch(@(dy) EntropyMetric(dy,kSpace,segkStart,segkStop,baseKSpaceRange(1),baseKSpaceRange(2),2),-a,options);
                % Find Entropy
                E = EntropyMetric(scout_dy,kSpace,segkStart,segkStop,baseKSpaceRange(1),baseKSpaceRange(2),2);
                if(isempty(minE)||E<minE)
                    % store dy and E
                    dy = scout_dy;
                    minE = E;
                end 
            end
        end
        % Shift y axis
        shift = exp(-i*2*pi/dims(1) .* [0:dims(1)-1] .* dy * dims(1))';
        shift = repmat(shift,[1 dims(1)]);
        shift = shift(segkStart:segkStop,:);
        kSpace(segkStart:segkStop,:) = kSpace(segkStart:segkStop,:).*shift;
        
        dx,dy
    end
    
    % Switch side of baseline
    positiveSide = mod(positiveSide+1,2);
    
    % if switch to positive side change base image
    if (positiveSide)
        baseKSpaceRange(1) = baseKSpaceRange(1)-2;
        baseKSpaceRange(2) = baseKSpaceRange(2)+2;
    end
    
end

corr_kSpace = kSpace;