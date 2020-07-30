function [D,J] = F_new(m,model,nw,flag)
% Forward operator
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

    %% size
    nr = length(model.zr);
    ns = length(model.zs);
    nx = prod(model.n);
    
	%% generate matrices
	A  = getA(model.f,m,model.h,model.n);
	P  = getP(model.h,model.n,model.zr,model.xr);
	if(flag == 1) % forward modeling
		Q  = getQ_for(model.h,model.n,model.zs,model.xs,nw,model.f);
	else          % inversoin
		Q  = getQ_inv(model.h,model.n,model.zs,model.xs,nw,model.f);
	end

	%% solve
	D = P'*(A\Q);
    D = D(:);
    
	%% Jacobian
	J = opFunction(nr*ns, nx, @(x,flag)Jmv(x,m,full(A\Q),model,flag));
	
end


function y = Jmv(x,m,U,model,flag)
    % size
    nr = length(model.zr);
    ns = length(model.zs);
    nx = prod(model.n);
    
    %% get matrices
    Pr = getP(model.h,model.n,model.zr,model.xr);
 
    %% compute mat-vec
    if flag == 1
         y = zeros(nr,ns);
        
         Rk = zeros(nx,ns);
         Ak = getA(model.f,m,model.h,model.n);
         Gk = @(u)getG(model.f,m,u,model.h,model.n);
         for l = 1:ns
            Rk(:,l) = -Gk(U(:,l))*x;
         end
         y(:,:) = Pr*(Ak\Rk);
      
        y = y(:);
    else
        y = zeros(nx,1);
        x = reshape(x,[nr,ns]);
       
        Ak = getA(model.f,m,model.h,model.n);
        Gk = @(u)getG(model.f,m,u,model.h,model.n);
        Rk = Ak'\(Pr*x(:,:));
        for l = 1:size(U,2)
            y = y - Gk(U(:,l))'*Rk(:,l);       
        end
    end

end
