function r0 = rtm(m0,model)
% rtm operator

%
%   D = P^TA^{-1}(m)Q
%
% where P, Q encode the receiver and source locations and L is the first-order FD matrix
%
% use:
%   [D,J] = F(m,model);
%
% input:
%   m - squared-slownes [s^2/km^2]
%   model.h - gridspacing in each direction d = [d1, d2];
%   model.n - number of gridpoints in each direction n = [n1, n2]
%   model.f - frequency [Hz].
%   model.{zr,xr} - {z,x} locations of receivers [m] (must coincide with gridpoints)
%   model.{zs,xs} - {z,x} locations of sources [m] (must coincide with gridpoints)
%
%
% output:
%   D - data matrix
%   J - Jacobian as Spot operator

%% generate matrices
P  = getP(model.h,model.n,model.zr,model.xr);
Q  = getQ_for(model.h,model.n,model.zs,model.xsf,model.nw,model.f);
  
%% rtm
r   = zeros(prod(model.n),1);
il  = zeros(prod(model.n),1);
nf  = 12;    f0 = 3;   
df  = 1;     k  = 1;

for f = f0:df:nf
     model.f = f;  k  = k+1;
     % real data
     tmp = ['D_' num2str(k) ];
     dobs = eval(tmp);
    
     % synthetic data
     A  = getA(model.f,m0,model.h,model.n);
     U  = A \ Q;
     D  = P' * U;  
        
     % RTM       
     V  = A' \ ( P * (D - dobs) );
     r  = r + real(sum(conj(U).*V,2));
     il = il + real(sum(conj(U).*U,2));
end

r0 = reshape(r,model.n) ./ reshape(il,model.n);
r0(1:2,:) = 0;
r0 = r0(:);
       
end

