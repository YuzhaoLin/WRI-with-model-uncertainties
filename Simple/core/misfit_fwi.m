function [f,g] = misfit_fwi(m,D,alpha,model)
%% Evaluate least-squares misfit
%
%   0.5||P^TA^{-1}(m)Q - D||_{F}^2 + 0.5\alpha||Lm||_2^2,
%
% where P, Q encode the receiver and source locations and L is the first-order FD matrix
%
% use:
%   [f,g,H] = misfit(m,D,alpha,model)
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
%   H - GN Hessian (function handle)


%% get matrices
m = m(:);
L = getL(model.h,model.n);
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
Q = getQ(model.h,model.n,model.zs,model.xs,model.f);
G = @(u)getG(model.f,u,model.n);

%% forward solve
U = A\Q;

%% compute f
f = .5*norm(P'*U - D,'fro')^2;

%% adjoint solve
V = A'\(P*(D - P'*U));

%% compute gradient
g = alpha*(L'*L)*m;

for k = 1:size(U,2)
    g = g + real(G(U(:,k))'*V(:,k));
end


end
