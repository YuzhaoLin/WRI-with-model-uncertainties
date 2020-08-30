
clc,clear
%% basic paramters
T  = 1;
dt = 1e-2;

c = 2000;

fs = 10;
ts = 0.2;

r = [800 1000 1200];

%% forward operator
params.T  = T;
params.dt = dt;
params.r  = r;

nr = length(params.r);
nt = (params.T/params.dt + 1);
F  = @(c)opFunction(length(params.r)*(params.T/params.dt + 1),(params.T/params.dt + 1),@(q,mode)Fmult(q,c,mode,params));

%% forward modeling
q = ricker(fs,ts,params);
d = F(c)*q;

%% calculate misfit function
rho = 1e-10;
cs = linspace(1500,2500,100);
phi = zeros(length(cs),2);
for k = 1:length(cs)
    phi(k,:) = misfit(cs(k),q,d,rho,params);
end

%% plot
figure;
plot(1e-3*cs,phi(:,1)/max(abs(phi(:,1))),'r','linewidth',2); hold on
plot(1e-3*cs,phi(:,2)/max(abs(phi(:,2))),'b','linewidth',2); 
xlabel('c [km/s]','fontsize',20);ylabel('\phi','fontsize',20);
set(gca,'fontsize',20)

%% save figures
print(1,'-depsc','wri-un-nonlinear.eps')
