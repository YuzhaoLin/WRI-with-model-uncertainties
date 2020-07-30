function [f,g] = misfit(m,D,alpha,model,Pp,Pm)
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
%   H - GN Hessian (function handle)


%% get matrices
mk = m(:);
L  = getL(model.h,model.n);
Ak = getA(model.f,mk,model.h,model.n);
P  = getP(model.h,model.n,model.zr,model.xr);
Q  = getQ_for(model.h,model.n,model.zs,model.xsf,model.nw,model.f);
G  = @(u)getG(model.f,mk,u,model.h,model.n);

%% forward solve
U0  = Ak \ Q;
D0  = P' * U0;

%% cpmpute PP
PM = P' * inv( Ak' * inv(Pp) * Ak ) * P;

%% compute f
f = .5*norm(D0 - D,'fro')^2 + .5*alpha*norm(L*mk)^2;

%% compute gradient residual 
K     = PM + Pm;
h0    = D - D0;
h_mdd =  K' * h0;
for iter = 1:3
    res = h0 - K*h_mdd;
    gh  = -K'*res;
    h_mdd = h_mdd - gh;
end
    
%% compute adjoint state
v0 = Ak' \ (P * h_mdd ); 
w0 = Ak \ (Pp * v0);

%% compute g
g1 = alpha*(L'*L)*mk;
g2 = alpha*(L'*L)*mk;

for k = 1:size(U0,2)
    g1 = g1 + 2*real(G(U0(:,k))'*v0(:,k)) ;
    g2 = g2 + 2*real(G(v0(:,k))'*w0(:,k)) ;
end
g = g1 - g2;
% tmp = reshape(g,model.n);
% tmp(1:5,:) = 0;
% g   = tmp(:);

%% get H
H = @(m)Hmv(m,U,alpha,model);
%H = H(m);

end

function y = Hmv(m,U,alpha,model)
%% get matrices
L = getL(model.h,model.n);
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% compute mat-vec
y = alpha*(L'*L)*m;

for k = 1:size(U,2);
    y = y + real(G(U(:,k))'*(A'\((P*P')*(A\(G(U(:,k))*m)))));
end

end
