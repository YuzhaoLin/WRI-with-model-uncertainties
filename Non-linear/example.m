
clc,clear
%% define model
n = [101 201];
h = [25 25];
[zz,xx] = ndgrid((0:n(1)-1)*h(1),(0:n(2)-1)*h(2));

v = @(v0,alpha)(v0 + alpha*zz(:));

%% make data
model.n = n;
model.h = h;
model.f = 5;
model.Is = {2,1};
model.Ir = {2,1:201};

% Noise 
sigm = 1e-20; 
mum    = zeros(1,length(model.Ir{2})); 
sigmm  = sigm*eye(length(model.Ir{2}),length(model.Ir{2}));
Mnoise = mvnrnd(mum, sigmm, length(model.Is{2}))';

Q = 1;
D = F(1./v(2000,0.75).^2,Q,model) + Mnoise;

%% scan over v0
vs = 1750:25:2250;

fv_fwi  = zeros(1,length(vs));
fv_new  = zeros(1,length(vs));

for k = 1:length(vs);
    vk = v(vs(k),0.75);
    
    fv_fwi(k)  = misfit_fwi(1./vk.^2,Q,D,model);
    fv_new(k)  = misfit_new(1./vk.^2,Q,D,model,sigmm);
  
end
   
%% plot results
figure;plot(vs,fv_fwi/max(fv_fwi),'linewidth',2); hold on
       plot(vs,fv_new/max(fv_new),'r','linewidth',2);
       xlabel('c [m/s]','fontsize',18);
       ylabel('Normalized misfit','fontsize',18);
       legend('FWI','New method','Location','SouthEast');
       set(gca,'fontsize',18); axis tight

%% save figures
print(1,'-depsc',['misfit']);
