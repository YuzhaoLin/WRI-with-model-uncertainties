%% setup


% read model
dx = 10;  n = [100, 100];
v  = 2*ones(n);              % initial model
m  = 1./(v(:)).^2;

% frequency
f = 5;   a = 1e-2;  b = 1e-3;

% receivers
xr = 20:1*dx:980;
zr = 2*dx*ones(1,length(xr));

% sources
xs = 20:10*dx:980;   
zs = 2*dx*ones(1,length(xs));

% regularization parameter
alpha = 0;   k = 1;

% grid
n  = size(v);          h  = dx*[1 1];
z  = [0:n(1)-1]*h(1);  x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

% parameters
model.n  = n;    model.h   = h;    
model.zr = zr;   model.xr  = xr;    model.f  = f;
model.zs = zs;   model.xs  = xs;    model.nf = 12;

%% Noise in source and data
sigm = 2e-10;  sigp = 1e-12;
% data covariance
mum    = zeros(1,length(xr)); 
sigmm  = sigm*eye(length(xr),length(xr));
Mnoise = mvnrnd(mum, sigmm, length(xs))';

% source covariance
mus    = zeros(1,n(1)*n(2)); 
sigmp  = sigp*eye(n(1)*n(2),n(1)*n(2));
Pnoise = mvnrnd(mus, sigmp, length(xs))';

%% operators
A  = magic(n(1)*n(2));
P  = getP(model.h,model.n,model.zr,model.xr);
Q  = getQ_for(model.h,model.n,model.zs,model.xs,model.nf,model.f);


%% random observed data
D = random('Normal',0,1,length(xr),length(xs));

%% FWI
tic;
U  =  A  \ Q;
V  =  A' \ (P*(D - P'*U));
toc;

%% WRI
tic;
U  = (a * A' * A + b * P * P') \ ( a * A' * Q  + b * P * D );
V  = (A * U - Q); 
toc;

%% New formula \sigma = I
tic;
U  =  A  \ Q;
H  =  P' * inv( A' * A ) * P;
V  =  A' \ (P*(D - P'*U));
W  =  A  \ V;
toc;

%% New formula \sigma = I
tic;
U  =  A  \ Q;
V  =  A' \ (P*(D - P'*U));
toc;



