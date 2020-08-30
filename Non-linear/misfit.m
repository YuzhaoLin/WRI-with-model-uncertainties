function f = misfit(c,q,d,rho,params)
%% calculation of misfit function for extended FWI

%% basic paramters
nr = length(params.r);
nt = (params.T/params.dt + 1);

%% forward operator
F = @(c)opFunction(length(params.r)*(params.T/params.dt + 1),(params.T/params.dt + 1),@(q,mode)Fmult(q,c,mode,params));

%% calculate misfit
dp = F(c)*q;
r  = d - dp;
K  = F(c)*F(c)';
s  = (rho*opDirac(nt*nr) + K)\r;
f  = [norm(d - dp)^2, r'*s];

end
