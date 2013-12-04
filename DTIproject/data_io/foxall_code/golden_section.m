%
%
% GOLDEN SECTION - Minimize a function in one dimension
%
% min_f   (string)     Name of a matlab function that returns values for the function to be minimized
% x_s     (float)      Start value of the interval in which the minimum value is expected to be located
% x_m     (float)      A point such that: min_f(x_s) > min_f(x_m) < min_f(x_e) and x_s < x_m < x_e
% x_e     (float)      End value of the interval in which the minimum is expected to be located.
% tol     (float)      A tolerence value for closing the interval
% p1 - p6 (float)      Parameter values to be passed to the function min_f
%
%

function [minx,minf] = golden_section(min_f,x_s,x_m,x_e,tol,p1,p2,p3,p4,p5,p6)

%----------------------------------------------------------------------
% CONSTANTS
%----------------------------------------------------------------------

R = 0.61803399;
C = 1-R;


%----------------------------------------------------------------------
% Function Evaluation Strings
%----------------------------------------------------------------------

h0 = [min_f,'(x0'];
h1 = [min_f,'(x1'];
h2 = [min_f,'(x2'];
h3 = [min_f,'(x3'];

if (nargin > 5)
    for I = 1:(nargin - 5)
       h0 = [h0,',p',int2str(I)];
       h1 = [h1,',p',int2str(I)];
       h2 = [h2,',p',int2str(I)];
       h3 = [h3,',p',int2str(I)];
    end
end

s0 = [h0,')'];
s1 = [h1,')'];
s2 = [h2,')'];
s3 = [h3,')'];



%----------------------------------------------------------------------
% INITIALIZATION
%----------------------------------------------------------------------

x0 = x_s;
x3 = x_e;
f0 = eval(s0);
f3 = eval(s3);

if (abs(x_e - x_m) > abs(x_m - x_s))
   x1 = x_m;
   x2 = x_m + C*(x_e - x_m);
else
   x2 = x_m;
   x1 = x_m - C*(x_m - x_s);
end

f1   = eval(s1);
f2   = eval(s2);
itry = 0;

%-----------------------------------------------------------------------
% DEBUG PRINT
%
fprintf(1,'%s\n','');
fprintf(1,'try = %d, f0 = %f, f1 = %f, f2 = %f, f3 = %f\n',itry,f0,f1,f2,f3);
fprintf(1,'try = %d, x0 = %f, x1 = %f, x2 = %f, x3 = %f\n',itry,x0,x1,x2,x3);
fprintf(1,'%s\n','');


%------------------------------------------------------------------------
% IMPROVE MINIMUM ESTIMATE
%------------------------------------------------------------------------

while( abs(x3-x0) > tol*(abs(x1)+abs(x2)) )

    if (f2 < f1)
       x0 = x1;
       x1 = x2;
       x2 = R*x1 + C*x3;
       f0 = f1;
       f1 = f2;
       f2 = eval(s2);
    else
       x3 = x2;
       x2 = x1;
       x1 = R*x2 + C*x0;
       f3 = f2;
       f2 = f1;
       f1 = eval(s1);
    end

itry = itry + 1;

fprintf(1,'try = %d, f0 = %f, f1 = %f, f2 = %f, f3 = %f\n',itry,f0,f1,f2,f3);
fprintf(1,'try = %d, x0 = %f, x1 = %f, x2 = %f, x3 = %f\n',itry,x0,x1,x2,x3);
fprintf(1,'%s\n','');

end

%--------------------------------------------------------------------------
% OUTPUT STAGE
%--------------------------------------------------------------------------

if (f1 < f2)
    minx = x1;
    minf = f1;
else
    minx = x2;
    minf = f2;
end

%-----------
% THE END
%-----------