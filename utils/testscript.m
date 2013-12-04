    
    % Question Set-up (current best
    xdim = 500; 
    ydim = 500; 
    n_timepoints = 10; % for example
    
    regressor = randn(xdim,ydim,n_timepoints); % 500x500x10
    % a stack of images w/ some interesting difference
    regressand = randn(xdim,ydim);  %500x500

    % Actual work  
    tic  % start timing
    % reshape to remove 2nd loop
    regressor = reshape(regressor, [xdim*ydim,n_timepoints])'; % not transpose
    regressand = regressand(:);
    % initialize
    est = zeros(size(regressor,2),1); 
    X = randn(n_timepoints,10);
    for x=1:size(picture,2) % second dimension b/c already transposed
        
        Y = regressand(x); % Nx1
        X(:,2)=picture(:,x); % Nx2
        
        XT = X'; % pre-transpose
        
        B = (XT*X)^(-1)*XT*Y; % B is 2x1 (always something >1 by 1)
        est(x) = B(1);
    end
    est = reshape(est,[xdim,ydim]); 
    toc % end timing

    
    
    % edited version
    % Question Set-up 
    tic
    xdim = 500; 
    ydim = 500; 
    n_timepoints = 10; % for example
    % Actual work
    picture = randn(xdim,ydim,n_timepoints);
    picture = reshape(picture, [xdim*ydim,n_timepoints])'; % note transpose
    est = zeros(size(picture,2),1); % initialization 
    for x=1:size(picture,2) % second dimension b/c already transposed
        X = picture(:,x); 
        Y = randn(n_timepoints,1);
        B = (X'*X)\X'*Y; 
        est(x) = B(1);

    end
    est = reshape(est,[xdim,ydim]); 
    toc

    
    
    % origianl
    tic
    xdim = 500; 
    ydim = 500; 
    n_timepoints = 10; % for example
    est = zeros(xdim,ydim,1); % initialization with explicit size
    picture = randn(xdim,ydim,n_timepoints);
    for x=1:size(picture,1)
        for y=1:size(picture,2)
            X = [squeeze(picture(x,y,:))]; % or similar creation of X matrix
            Y = randn(n_timepoints,1);
            B = (X'*X)^-1*X'*Y; % SEE HERE 
            % this will not be a scalar. It will be size(X,2) by 1
            % sometimes you keep everything and do
            % estimate(x,y,:) = B(:);
            % sometimes just the first element is important and you do
            est(x,y) = B(1);
        end
    end
    toc
    
    
    
    % using permmuite
    tic
    xdim = 500; 
    ydim = 500; 
    n_timepoints = 10; % for example
    est = zeros(xdim,ydim,1); % initialization with explicit size
    picture = randn(xdim,ydim,n_timepoints);
    picture = permute(picture,[3 1 2]);
    for x=1:size(picture,2)
        for y=1:size(picture,3)
            X = picture(:,x,y); % or similar creation of X matrix
            Y = randn(n_timepoints,1);
            B = (X'*X)^-1*X'*Y; % SEE HERE 
            % this will not be a scalar. It will be size(X,2) by 1
            % sometimes you keep everything and do
            % estimate(x,y,:) = B(:);
            % sometimes just the first element is important and you do
            est(x,y) = B(1);
        end
    end
    toc