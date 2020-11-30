function P = getQ(h,n,zt,xt,iw)
%% Define sampling operator
%
% use:
%   P = getP(h,n,zt,xt,nf,iw)
%
% input:
%   h,n   - gridspacing and number of gridpoints
%   zt,xt - arrays defining sampling points (must coincide with grid)
%   iw - frequency infor to generate source wavelet 
%
% output
%   P     - sparse matrix

%% Ricker wavelet
fm = 8;   
dt = 0.001;
t  = 0:1:499;  

tmp    = pi*fm*(t-200)*dt;
ricker = (1-2*(tmp).^2).*exp(-(tmp).^2);

% frequency domain
rf = fft(ricker);

z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

for k = 1:length(zt)
    i(k) = find((zz(:)==zt(k))&(xx(:)==xt(k)));
end

I = speye(prod(n));
P = abs(rf(floor(iw)))*I(:,i)/sqrt(prod(h));



end