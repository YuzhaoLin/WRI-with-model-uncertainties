function [f,g] = misfit_fwiqq(m,D,alpha,model,Pm)
%% Evaluate least-squares misfit
%
%   0.5||P^TA^{-1}(m)Q - D||_{F}^2_{Sigma(m)}
%
% where P, Q encode the receiver and source locations and L is the first-order FD matrix
%
% use:
%   [f,g] = misfit(m,D,model)
%
% input:
%   m - squared-slownes [s^2/km^2]
%   D - single-frequency data matrix
%   alpha - regularization parameter
%   model.h - gridspacing in each direction d = [d1, d2];
%   model.n - number of gridpoints in each direction n = [n1, n2]
%   model.f - frequency [Hz].
%   model.{zr,xr} - {z,x} locations of receivers [m] (must coincide with gridpoints)
%   model.{zs,xs} - {z,x} locations of sources [m] (must coincide with gridpoints)
%
%
% output:
%   f - value of misfit
%   g - gradient (vector of size size(m))


%% get matrices
mk = m(:);
L  = getL(model.h,model.n);
Ak = getA(model.f,mk,model.h,model.n);
P  = getP(model.h,model.n,model.zr,model.xr);
Q  = getQ_for(model.h,model.n,model.zs,model.xs,model.nf,model.f);
G  = @(u)getG(model.f,u,model.n);

%% forward solve
U0  = Ak \ Q;
D0  = P' * U0;

%% compute weighted residual by the Sherman¨CMorrison formula
PM    = Pm; 
weigi = zeros(length(model.xr),length(model.xr));
for is = 1:size(U0,2)
    weigi = weigi + PM - (PM * D0(:,is) * D0(:,is)' * PM)/(1 + D0(:,is)' * PM * D0(:,is)) ;
end
h_sm = weigi * (D-D0);

%% compute f
f  = .5*norm( h_sm .* (D-D0) )^2;

%% compute adjoint state
v0  = Ak' \ (P  * h_sm ); 
w0  = Ak \ Q * Q' * v0;

%% compute g
g1 = alpha*(L'*L)*mk;
g2 = alpha*(L'*L)*mk;

for k = 1:size(U0,2)
    g1 = g1 + 2*real(G(U0(:,k))'*v0(:,k)) ;
    g2 = g2 + 2*real(G(v0(:,k))'*w0(:,k)) ;    
end
g  = g1 - g2;


end
