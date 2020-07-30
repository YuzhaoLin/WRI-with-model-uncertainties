function P = getQ_inv(h,n,zt,xt,nw,iw)
% Define sampling operator
%
% use:
%   P = getP(h,n,zt,xt)
%
% input:
%   h,n   - gridspacing and number of gridpoints
%   zt,xt - arrays defining sampling points (must coincide with grid)
%
% output
%   P     - sparse matrix

% wavelet
f = 8; dt = 0.01;  nc = floor(nw/2);
alpha = -4*f*f*pi*pi/log(0.1);
beta  = (nc-iw+5);
w     = -8*pi*pi*beta*dt*exp(-alpha*dt*dt*beta*beta);

% coordinates
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

for k = 1:length(zt)
    i(k) = find((zz(:)==zt(k))&(xx(:)==xt(k)));
end

I = speye(prod(n));
P = abs(w)*I(:,i)/sqrt(prod(h));
