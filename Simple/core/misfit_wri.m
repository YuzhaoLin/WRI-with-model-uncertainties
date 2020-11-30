function [f,g] = misfit_wri(m,D,alpha,model,SP,SM)
%% Evaluate least-squares misfit
%
%   0.5||P^U - D||^2_{\Sigma_m} + 0.5||AU-Q||^2_{\Sigma_p},
%
% where P, Q encode the receiver and source locations and L is the first-order FD matrix
%
% use:
%   [f,g] = misfit(m,D,alpha,model,SP,SM)
%
% input:
%   m - squared-slownes [s^2/km^2]
%   D - single-frequency data matrix
%   alpha - regularization parameter
%   model.h - gridspacing in each direction d = [d1, d2];
%   model.n - number of gridpoints in each direction n = [n1, n2]
%   model.f - frequency [Hz].
%   SP,SM   - theory and measurement covariance matrix
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
Q  = getQ(model.h,model.n,model.zs,model.xs,model.f);
G  = @(u)getG(model.f,u,model.n);

%% wavefield reconstruction
U  = (Ak' *  Ak / SP + P * P' / SM) \ ( Ak' * Q / SP  + P * D  / SM);

%% compute adjoint state
v0 =  (Ak * U - Q); 

%% compute f
f = .5*norm(P' * U - D,'fro')^2 + .5*norm(v0,'fro')^2;

%% compute g
g  = alpha*(L'*L)*mk;

for k = 1:size(U,2)
     g  = g + real(G(U(:,k))'*v0(:,k)) ;
end


end

