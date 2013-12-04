function c=vuCubicBSplineLineFit(y)
% Cubic B-Spline line fitting
% Adapted for vuTools : April 7, 2008 by: Kevin Wilson
% Based on work by Biomedical Imaging Group at EPA
%
% Given values (y) of a function at x=1,2,...n calculate coefficients (c)
% of a B-spline interpolation
%
% INPUT :
%   y : Value of function at 1,2,..n
%
% Output :
%   c : Cubic B-Spline coeffiecients
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

if(nargin<1)
    error('MATLAB:vuCubicBSplineLineFit:NotEnoughInputs', 'Not enough input arguments.');
end

N=length(y);
A=zeros(N);
y=y(:);

if N==1
    c=1.5*y;
else
    
  A(1,1:2)=[2/3 1/6];
  A(N,N-1:N)=[1/6 2/3];
  
  for i=2:N-1
    A(i,i-1:i+1)=[1/6 2/3 1/6];
  end  

  c=A\y;
  c=c(:)';
  
end 