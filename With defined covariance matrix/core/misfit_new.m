function [f,g,H] = misfit_new(m,D,alpha,model,nw)
% Evaluate least-squares misfit
%
%   0.5||P^TA^{-1}(m)Q - D||_{F}^2 + 0.5\alpha||Lm||_2^2,
%
% where P, Q encode the receiver and source locations and L is the first-order FD matrix
%
% use:
%   [f,g,H] = misfit(m,D,model)
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
%   H - GN Hessian (Spot operator)


%% get matrices
m  = m(:);
L  = getL(model.h,model.n);
Ap = getA(model.f,m,model.h,model.n);
P  = getP(model.h,model.n,model.zr,model.xr);
Q  = getQ_for(model.h,model.n,model.zs,model.xs,nw,model.f);

%% forward modeling
Dp = P'*(Ap\Q);
Dp = Dp(:);
    
%% Jacobian
Jp = opFunction(length(model.zr)*length(model.zs), prod(model.n), @(x,flag)Jmv(x,m,full(A\Q),model,flag));

%% compute covariance 
sp = gama*(Q')*(Q')';
sm = (D-Dp) * (D-Dp)';

%% compute f
f = .5*norm(Dp - D,'fro')^2 + .5*alpha*norm(L*m)^2;

%% compute r
r   = ( (Ap'\P * sm * P*Ap') + sp ) \ ( Ap'\P * sm * (D-Dp) );

%% compute gradient g
u0 = Ap\Q;
w0 = gama*(Q')*(Q')' * (Ap')*(Ap')*r;
v0 = (P*Ap\)' * (sm + sp) * (D-Dp);
g1 = -2*u0' * (Jp') * v0;
g2 = 2*v0' * (Jp') * w0;
g  = g1 + g2;

%% get H
H = Jp'*Jp + alpha*opMatrix(L'*L);

end
