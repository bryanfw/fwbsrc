%im = IMG;
out = im;
K = 10;
dt = 0.125;

c1 = zeros(size(im,1),size(im,2));
c2 = zeros(size(im,1),size(im,2));
c3 = zeros(size(im,1),size(im,2));
c4 = zeros(size(im,1),size(im,2));
c5 = zeros(size(im,1),size(im,2));
c6 = zeros(size(im,1),size(im,2));
c7 = zeros(size(im,1),size(im,2));
c8 = zeros(size(im,1),size(im,2));

for iter = 1:10
    
    % Loop through each pixel
    for i = 1:size(im,1)
        for j = 1:size(im,2)
            
            % Reset Derivatives
            d1 = 0;
            d2 = 0;
            d3 = 0;
            d4 = 0;
            d5 = 0;
            d6 = 0;
            d7 = 0;
            d8 = 0;

            for k = 1:size(im,3)
                
                % Caluculate Derivatives in each direction
                try
                    d1 = d1 + (out(i,j,k) - out(i-1,j-1,k))^2;
                catch
                    d1 = 0;
                end
                try
                    d2 = d2 + (out(i,j,k) - out(i-1,j,k))^2;
                catch
                    d2 = 0;
                end
                try
                    d3 = d3 + (out(i,j,k) - out(i-1,j+1,k))^2;
                catch
                    d3 = 0;
                end
                try
                    d4 = d4 + (out(i,j+1,k) - out(i,j,k))^2;
                catch
                    d4 = 0;
                end
                try
                    d5 = d5 + (out(i+1,j+1,k) - out(i,j,k))^2;
                catch
                    d5 = 0;
                end
                try
                    d6 = d6 + (out(i+1,j,k) - out(i,j,k))^2;
                catch
                    d6 = 0;
                end
                try
                    d7 = d7 + (out(i+1,j-1,k) - out(i,j,k))^2;
                catch
                    d7 = 0;
                end
                try
                    d8 = d8 + (out(i,j,k) - out(i,j-1,k))^2;
                catch
                    d8 = 0;
                end
            end
            
            % Calculate Curvatures
            c1(i,j) = exp(-1*((sqrt(d1)/K)^2));
            c2(i,j) = exp(-1*((sqrt(d2)/K)^2));
            c3(i,j) = exp(-1*((sqrt(d3)/K)^2));
            c4(i,j) = exp(-1*((sqrt(d4)/K)^2));
            c5(i,j) = exp(-1*((sqrt(d5)/K)^2));
            c6(i,j) = exp(-1*((sqrt(d6)/K)^2));
            c7(i,j) = exp(-1*((sqrt(d7)/K)^2));
            c8(i,j) = exp(-1*((sqrt(d8)/K)^2));

        end
    end

    
    for i = 1:size(im,1)
        for j = 1:size(im,2)
            for k = 1:size(im,3)
                
                % Caluculate Derivatives in each direction
                try
                    d1 = out(i,j,k) - out(i-1,j-1,k);
                catch
                    d1 = 0;
                end
                try
                    d2 = out(i,j,k) - out(i-1,j,k);
                catch
                    d2 = 0;
                end
                try
                    d3 = out(i,j,k) - out(i-1,j+1,k);
                catch
                    d3 = 0;
                end
                try
                    d4 = out(i,j+1,k) - out(i,j,k);
                catch
                    d4 = 0;
                end
                try
                    d5 = out(i+1,j+1,k) - out(i,j,k);
                catch
                    d5 = 0;
                end
                try
                    d6 = out(i+1,j,k) - out(i,j,k);
                catch
                    d6 = 0;
                end
                try
                    d7 = out(i+1,j-1,k) - out(i,j,k);
                catch
                    d7 = 0;
                end
                try
                    d8 = out(i,j,k) - out(i,j-1,k);
                catch
                    d8 = 0;
                end
                
                % Calculate Flows
                f1 = c1(i,j)*d1;
                f2 = c2(i,j)*d2;
                f3 = c3(i,j)*d3;
                f4 = c4(i,j)*d4;
                f5 = c5(i,j)*d5;
                f6 = c6(i,j)*d6;
                f7 = c7(i,j)*d7;
                f8 = c8(i,j)*d8;

                % Delta
                delta = 0.5*(f5 - f1) + f6 - f2 + 0.5*(f7 - f3) + f4 - f8;

                % Calculate new pixel value
                out(i,j,k) = out(i,j,k) + dt*delta;
            end
        end
    end

end