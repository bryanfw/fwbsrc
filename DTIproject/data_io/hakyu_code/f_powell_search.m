
function xi = f_powell_search(optfun,p,varargin)
%[f_powell_search] searches using Powell's method.
%
% Usage:
%   xi = f_powell_search(optfun,p,varargin)
%
% Input:
%   optfun -
%       Cost function for optimization
%   p -
%       Searching parameter
%   varargin -
%       {1} img_ref, reference image, to
%       {2} img_tmp, temporary image, from
%       {3} method_s, registration method, 'rigid' or 'affine'
%       {4} cost_s, cost function, 'nmi','mi','cor'
%           'nmi'   normalized mutual information
%           'mi'    mutual information
%           'ecc'   entropy correlation coefficient
%           'ncc'   normalized cross correlation
%           'cor'   correlation coefficient
%       {5} nBin, number of bins used for 'nmi' or 'mi'
%       {6} flag_display, flag for displaying parameter, 
%           'on' or 'off'    
%
% Output:
%   xi -
%       Final parameter searched
%
% Note:
%   This routine is from "Numerical Recipes in C" chap 10 and from
%   [spm_mireg.m].
%
%
% Last modified
% 2011.01.24.
%   Incorporated old scripts and functions in
%   'N:\Computer_DTI\C\home\Course\MATH287_05SP\Project'
%   Take 'flag_display' as 8th input.
% 2011.01.25.
%   Generate 'New' procedure.
% HKJ



%% Check input
if nargin == 6
    img_ref = varargin{1};
    img_tmp = varargin{2};
    method_s = varargin{3};
    cost_s = varargin{4};
    nBin = 64; % default
    flag_display = 'off'; % default
elseif nargin == 7
    img_ref = varargin{1};
    img_tmp = varargin{2};
    method_s = varargin{3};
    cost_s = varargin{4};
    nBin = varargin{5};
    flag_display = 'off'; % default
elseif nargin == 8
    img_ref = varargin{1};
    img_tmp = varargin{2};
    method_s = varargin{3};
    cost_s = varargin{4};
    nBin = varargin{5};
    flag_display = varargin{6};
else
    error('f_powell_search:main','Number of input must be 6 or 7')
end

% Check dimension.
DIM_ref = length(size(img_ref));
DIM_tmp = length(size(img_tmp));
if DIM_ref~=DIM_tmp
    error('f_powell_search:main','Dimension of the two images must be the same')
else
    DIM = DIM_ref; % or DIM_tmp
    if DIM > 3
        error('f_powell_search:main','Dimension of the two images must be 2-D or 3-D')
    end
end

% Check size.
if any(size(img_ref)~=size(img_tmp))
    error('f_powell_search:main','Size of the two images must be the same')
else
    [ny,nx,nz] = size(img_ref);
end


% Check search variables.
if size(p,2) ~= 1
    p = p';
end
xl = eye(length(p));    % line search direction
sc = ones(length(p),1); % parameter multiplier



% Run search.

%---------------------
% Old.
% for samp=1:size(sc,1),
%     [xi,fval,xl] = powell(p(:),xl,1e-3,optfun,img_ref,img_tmp,method_s,cost_s,nBin,flag_display);
%     xi = (xi.*sc)';
% end;


%---------------------
% New.

% Smoothing kernel.
kern_size_min = [12,12]; % min kernel size
kern_size_max = [24,24]; % max kernel size
kern_size = [nBin/2,nBin/2];
if any(kern_size <= kern_size_min)
    kern_size = kern_size_min;
else
    kern_size = kern_size_max;
end
fwhm = 7; % isotropic
smooth_kern = f_gauss_kernel(kern_size,fwhm); % 2-D kernel

% Re-sampling factor.
samp_factor = [0.5,1]; % over- (>1) and under- (<1) sampling factor

% Search.
for samp=1:length(samp_factor)
    
    % Report.
    if strcmpi(flag_display,'on')
        fprintf('\n\n')
        fprintf('Sampling factor [%f]\n',samp_factor(samp))
        fprintf('\n')
    end
    
    % Smooth and re-sample.
    if DIM==2
        img_ref1 = imresize(conv2(img_ref,smooth_kern,'same'),samp_factor(samp));
        img_tmp1 = imresize(conv2(img_tmp,smooth_kern,'same'),samp_factor(samp));
    else
        img_ref1 = zeros(ceil(size(img_ref)*samp_factor(samp)),'single');
        img_tmp1 = zeros(ceil(size(img_tmp)*samp_factor(samp)),'single');
        for indz = 1:nz
            img_ref1(:,:,indz) = imresize(conv2(img_ref(:,:,indz),smooth_kern,'same'),samp_factor(samp));
            img_tmp1(:,:,indz) = imresize(conv2(img_tmp(:,:,indz),smooth_kern,'same'),samp_factor(samp));
        end
    end
    
    [xi,fval,xl] = powell(p,xl,1e-3,optfun,img_ref1,img_tmp1,method_s,cost_s,nBin,flag_display);
    p = xi;
    %xi = (xi.*sc)';
end





%% Subfunctions

%__________________________________________________________________________
function [p,fret,xl] = powell(p,xl,ftol,func,varargin)
% Powell optimisation method - taken from Numerical Recipes (p. 417) and
% modified slightly.

p=p(:);
ITMAX = 32;
fret  = feval(func,p,varargin{:});
pt    = p;
for iter=1:ITMAX,
	fp   = fret;
	ibig = 0;
	del  = 0.0;
	for i=1:length(p),
		fptt = fret;
		[p,xit,fret] = linmin(p,xl(:,i),func,varargin{:});
		if abs(fptt-fret) > del,
			del  = abs(fptt-fret);
			ibig = i;
		end;
	end;
	if 2.0*abs(fp-fret) <= ftol*(abs(fp)+abs(fret)),
		return;
	end;
	ptt  = 2.0*p-pt;
	xit  = p-pt;
	pt   = p;
	fptt = feval(func,ptt,varargin{:});
	if fptt < fp,
		t = 2.0*(fp-2.0*fret+fptt)*(fp-fret-del).^2-del*(fp-fptt).^2;
		if t < 0.0,
			[p,xit,fret] = linmin(p,xit,func,varargin{:});
			xl(:,ibig)   = xl(:,end);
			xl(:,end)    = xit;
		end;
	end;
end;
warning('Too many iterations in routine POWELL');
return;





%__________________________________________________________________________
function [p,xl,fret] = linmin(p,xl,func,varargin)
% Code based on Numerical Recipes in C (p. 419)

global lnm
lnm = struct('pcom',p,'xicom',xl,'func',func,'args',[]);
lnm.args = varargin;
ax    = 0.0;
xx    = 1.0;
bx    = 2.0;
[ax,xx,bx,fa,fx,fb] = mnbrak(ax,xx);
[fret,xmin] = brent(ax,xx,bx,fx,2.0e-3);
xl    = xl * xmin;
p     = p + xl;
return;





%__________________________________________________________________________
function f = f1dim(x)
% Code based on Numerical Recipes in C (p. 419)

global lnm
xt = lnm.pcom+x.*lnm.xicom;
f = feval(lnm.func,xt,lnm.args{:});
return;





%__________________________________________________________________________
function [ax,bx,cx,fa,fb,fc] = mnbrak(ax,bx)
% Code based on Numerical Recipes in C (p. 400)
GOLD   = 1.618034;
GLIMIT = 100.0;
TINY   = 1.0e-20;

fa=f1dim(ax);
fb=f1dim(bx);

if fb > fa
	dum = ax; ax = bx; bx = dum;
	dum = fb; fb = fa; fa = dum;
end;
cx = bx+GOLD*(bx-ax);
fc = f1dim(cx);
while fb > fc,
	r    = (bx-ax)*(fb-fc);
	q    = (bx-cx)*(fb-fa);
	u    = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*(abs(q-r)+TINY)*sign(q-r));
	ulim = bx+GLIMIT*(cx-bx);
	if (bx-u)*(u-cx) > 0.0,
		fu=f1dim(u);
		if fu < fc,
			ax = bx; bx =  u;
			fa = fb; fb = fu;
			return;
		elseif fu > fb,
			cx = u;
			fc = fu;
			return;
		end;
		u  = cx+GOLD*(cx-bx);
		fu = f1dim(u);
	elseif (cx-u)*(u-ulim) > 0.0
		fu=f1dim(u);
		if fu < fc,
			bx = cx; cx = u; u = cx+GOLD*(cx-bx);
			fb = fc; fc = fu; fu = f1dim(u);
		end;
	elseif (u-ulim)*(ulim-cx) >= 0.0,
		u  = ulim;
		fu = f1dim(u);
	else,
		u  = cx+GOLD*(cx-bx);
		fu = f1dim(u);
	end;
	ax = bx; bx = cx; cx = u;
	fa = fb; fb = fc; fc = fu;
end;
return;





%__________________________________________________________________________
function [fx, x] = brent(ax,bx,cx,fx, tol)
% Code based on Numerical Recipes in C (p. 404)

ITMAX = 100;
CGOLD = 0.3819660; % 1-(1-sqrt(5))/2
e = 0.0;
a = min(ax,cx);
b = max(ax,cx);
x = bx; w = bx; v = bx;
fw = fx;
fv = fx;
for iter=1:ITMAX,
	xm   = 0.5*(a+b);
	tol1 = 2e-4*abs(x)+eps;
	tol2 = 2.0*tol1;
	if abs(x-xm) <= tol,
		return;
	end;
	if abs(e) > tol1,
		r     = (x-w)*(fx-fv);
		q     = (x-v)*(fx-fw);
		p     = (x-v)*q-(x-w)*r;
		q     = 2.0*(q-r);
		if q > 0.0, p = -p; end;
		q     = abs(q);
		etemp = e;
		e     = d;
		if abs(p) >= abs(0.5*q*etemp) | p <= q*(a-x) | p >= q*(b-x),
			if x >= xm, e = a-x; else, e = b-x; end;
			d = CGOLD*(e);
		else,
			d = p/q;
			u = x+d;
			if u-a < tol2 | b-u < tol2,
				d = tol1*sign(xm-x);
			end;
		end;
	else,
		if x>=xm, e = a-x; else, e = b-x; end;
		d = CGOLD*e;
	end;
	if abs(d) >= tol1, u = x+d; else, u = x+tol1*sign(d); end;
	fu=f1dim(u);
	if fu <= fx,
		if u >= x, a=x; else, b=x; end;
		 v =  w;  w =  x;  x =  u;
		fv = fw; fw = fx; fx = fu;
	else,
		if u < x, a=u; else, b=u; end;
		if fu <= fw | w == x,
			 v  = w;  w =  u;
			fv = fw; fw = fu;
		elseif fu <= fv | v == x | v == w,
			 v =  u;
			fv = fu;
		end;
	end;
end;
warning('Too many iterations in BRENT');
return;



