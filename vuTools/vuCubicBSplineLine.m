function yy=vuCubicBSplineLine(c,xx)
% Cubic B-Spline line
% Adapted for vuTools : April 7, 2008 by: Kevin Wilson
% Based on work by Biomedical Imaging Group at EPA
%
% Given coefficients (c) of a cubic B-spline calculate values (y) at
% locations (x)
%
% INPUT :
%   c : Cubic B-Spline coeffiecients
%   xx : Values to be calcuated at
%
% Output :
%   yy : Cubic B-Spline line value
%
% Example
%   x = 0:10;
%   y = sin(x);
%   xx = 1:1:10;
%   yy = spline(x,y,xx);
%   c = vuCubicBSplineLineFit(yy);
%   xNew = 0:.25:10;
%   yNew = vuCubicBSplineLine(c,xNew);
%   plot(x,y,'o',xNew,yNew)

if(nargin<2)
    error('MATLAB:vuCubicBSplineLine:NotEnoughInputs', 'Not enough input arguments.');
end

yy = zeros(1,length(xx));

for i = 1:length(xx)
    x = xx(i);

    lenc=length(c);
    xf=floor(x)-1; 

    if xf>=1 && xf<=lenc

      y=c(xf)*cubicBSpline(x-xf);

    else 

      y=0;

    end

    if xf+1>=1 && xf+1<=lenc

      y=y+c(xf+1)*cubicBSpline(x-xf-1);

    end 

    if xf+2>=1 && xf+2<=lenc

      y=y+c(xf+2)*cubicBSpline(x-xf-2);

    end 

    if xf+3>=1 && xf+3<=lenc

      y=y+c(xf+3)*cubicBSpline(x-xf-3);

    end
    
    yy(i) = y;
    
end


function y=cubicBSpline(x)
% Calculate the value of a cubic B-spline at point x

x=abs(x);

if x>2
    
  y=0;
  
else
    
  if x>1
      
    y=(2-x)^3/6;
    
  else
      
    y=2/3-x^2*(1-x/2);
    
  end
  
end