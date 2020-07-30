function P = getQ_for(h,n,zt,xt,nw,iw)
%% Define sampling operator
%
% use:
%   P = getP(h,n,zt,xt,nw,iw)
%
% input:
%   h,n   - gridspacing and number of gridpoints
%   zt,xt - arrays defining sampling points (must coincide with grid)
%   nw,iw - frequency infor to general wavelet
%
% output
%   P     - sparse matrix

%% wavelet
f = 8; dt = 0.01;  nc = floor(nw/2); % f, dt can change based on wavelet
alpha = (nc-iw+1)*f*dt*pi;
beta  = alpha^2;
w     = (1-beta*2)*exp(-beta);

%% coordinates
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

for k = 1:length(zt)
    i(k) = find((zz(:)==zt(k))&(xx(:)==xt(k)));
end

I = speye(prod(n));
P = abs(w)*I(:,i)/sqrt(prod(h));
