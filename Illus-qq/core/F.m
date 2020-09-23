function D = F(m,model,Pnoise)
%% Forward operator
%
%   D = P^TA^{-1}(m)(Q+Pnoise)
%
% where P, Q encode the receiver and source locations and L is the first-order FD matrix
%
% use:
%   D = F(m,model,nw,Pnoise);
%
% input:
%   m - squared-slownes [s^2/km^2]
%   model.h - gridspacing in each direction d = [d1, d2];
%   model.n - number of gridpoints in each direction n = [n1, n2]
%   model.f - frequency [Hz] model.nw - max frequency
%   model.{zr,xr} - {z,x} locations of receivers [m] (must coincide with gridpoints)
%   model.{zs,xs} - {z,x} locations of sources [m] (must coincide with gridpoints)
%   Pnoise - source noise represnets foward model uncertainties
%
%
% output:
%   D - data matrix

%% generate matrices
A  = getA(model.f,m,model.h,model.n);
P  = getP(model.h,model.n,model.zr,model.xr);
Q  = getQ_for(model.h,model.n,model.zs,model.xsf,model.nf,model.f);

%% solve
D = P'* (A\(Q+Pnoise));
