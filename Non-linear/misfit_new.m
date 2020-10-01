function [f] = misfit_pen(mt,Q,D,model,Pm)
%% Penalty objective:
%	
%	.5*||Pu - d||^2 + .5\lambda^2||A(m)u - q||^2,
% 
%  where 
%   
%   u = (P^*P + \lambda^2A^*A)^{-1}(P^*d + \lambda^2A^*q). 
%
%% use:
%    [f,g,h] = misfit_pen(mt,Q,D,lambda,model)

%% input:
%   mt - model [s^2/m^2]
%   Q  - sources
%   D  - data
%   lambda - penalty parameter
%   model.h - [dz,dx] gridspacing in z and x direction [m]
%   model.n - [nz,nx] number of gridpoints in z and x direction
%   model.f - frequencies
%   model.Is - cell array with source location indices {Iz, Ix}
%   model.Ir - cell array with receiver location indices {Iz,Ix}

%% output:
%	f - value
%   g - gradient
%   h - diagonal of Hessian

%% This program is part of the paper
% "Mitigating local minima in full-waveform inversion by expanding the search space",
% T. van Leeuwen and F.J. Herrmann, 2013 (submitted to GJI).
%
% Copyright (C) 2013 Tristan van Leeuwen (tleeuwen@eos.ubc.ca)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% get sampling operators
Ps = getP(model.n,model.Is{:});
Pr = getP(model.n,model.Ir{:});

%% initialize misfit, gradient and Hessian
f = 0;

%% loop over frequencies
for k = 1:length(model.f) 
	% get Helmholtz operator
	At = getA(model.f(k),mt(:),model.h,model.n);

	% solve wave-equation
	Ut = At\Ps'*Q;
	Dt = Pr*Ut;
    
	% compute weighted residual by the Sherman¨CMorrison formula
    PM    = Pm; 
    weigi = zeros(length(model.Ir{2}),length(model.Ir{2}));
    for is = 1:size(Ut,2)
        weigi = weigi + PM - (PM * Dt(:,is) * Dt(:,is)' * PM)/(1 + Dt(:,is)' * PM * Dt(:,is)) ;
    end
    h_sm = weigi * (D-Dt);

    % compute f
    f  = .5*norm( h_sm .* (D-Dt) )^2;
end