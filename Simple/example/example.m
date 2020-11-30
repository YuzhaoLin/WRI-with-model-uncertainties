%% setup


% read model
dx = 10;  n = [81, 121];
v0 = 2*ones(n);              % initial model
dv = zeros(n);
dv(33:47,53:67) = 0.2;       % abnormal
dv(66:69,:)     = 0.3;     

m  = 1./(v0(:) + dv(:)).^2;
mk = 1./(v0(:)).^2;

% set frequency
nf = 12;    f0 = 3;   df = 1;

% receivers
xr = 20:1*dx:1190;
zr = 2*dx*ones(1,length(xr));

% sources
xs = 20:10*dx:1190;   
zs = 2*dx*ones(1,length(xs));

% regularization parameter
alpha = 0;   k = 1;

% grid
n  = size(v0);         h  = dx*[1 1];
z  = [0:n(1)-1]*h(1);  x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

% parameters
model.n  = n;    model.h   = h;    
model.zr = zr;   model.xr  = xr;    model.nf  = nf;
model.zs = zs;   model.xs  = xs;    model.dx  = dx;

%% Noise in source and data
sigm = 5e-8;  sigp = 5e-8;
% data covariance
mum    = zeros(1,length(xr)); 
sigmm  = sigm*eye(length(xr),length(xr));
Mnoise = mvnrnd(mum, sigmm, length(xs))';

% source covariance
mus    = zeros(1,n(1)*n(2)); 
sigmp  = sigp*eye(n(1)*n(2),n(1)*n(2));
Pnoise = mvnrnd(mus, sigmp, length(xs))';

%% forward
for f = f0:df:nf
    model.f = f;
    fprintf('Forward frequency: %3d \n', f);
    d      = F(m,model,Pnoise); 
    d      = d + Mnoise;
    tmp    = ['D_' num2str(k) ];
    eval([tmp,'=d;']);  k = k+1;
end

%% inversion
k = 1; fwimk = 1./(v0(:)).^2;
for f = f0:df:nf
    tic;
    model.f = f;   
    fprintf('FWI Inversion frequency: %3d \n', f);
    % misfit
    tmp = ['D_' num2str(k) ];
    dobs = eval(tmp);
    fh = @(m)misfit_fwi(m,dobs,alpha,model);             % FWI
    %fh = @(m)misfit_wri(m,dobs,alpha,model,sigp,sigm);   % WRI
    %fh = @(m)misfit_fwii(m,dobs,alpha,model,sigp,sigmm); % FWI with Identity covariance
    %fh = @(m)misfit_fwiqq(m,dobs,alpha,model,sigmm);     % FWI with qq^* covariance
    %fh = @(m)misfit_fwiai(m,dobs,alpha,model,sigmm);     % FWI with source distance annihilator
   
    % Simple BB iteration
    [fwimk,hist] = BBiter(fh,fwimk,1e-20,5);   
	dlmwrite(['../input/hist_fwi_' num2str(f) '.txt'],hist);	
 
    k = k+1;   toc;    
end
vfwi = reshape(real(1./sqrt(fwimk)),n);

k = 1; wrimk = 1./(v0(:)).^2;
for f = f0:df:nf
    tic;
    model.f = f;   
    fprintf('FWI Inversion frequency: %3d \n', f);
    % misfit
    tmp = ['D_' num2str(k) ];
    dobs = eval(tmp);
    fh = @(m)misfit_wri(m,dobs,alpha,model,sigp,sigm);   % WRI
    
    % Simple BB iteration
    [wrimk,hist] = BBiter(fh,wrimk,1e-20,5);   
	dlmwrite(['../input/hist_wri_' num2str(f) '.txt'],hist);	
 
    k = k+1;   toc;    
end
vwri = reshape(real(1./sqrt(wrimk)),n);

% k = 1; fwiai = 1./(v0(:)).^2;
% for f = f0:df:nf
%     tic;
%     model.f = f;   
%     fprintf('FWI Inversion frequency: %3d \n', f);
%     % misfit
%     tmp = ['D_' num2str(k) ];
%     dobs = eval(tmp);
%     fh = @(m)misfit_fwiai(m,dobs,alpha,model,sigm);     % FWI with source distance annihilator
%    
%     % Simple BB iteration
%     [fwiai,hist] = BBiter(fh,fwiai,1e-20,5);   
% 	dlmwrite(['../input/hist_fwi_' num2str(f) '.txt'],hist);	
%  
%     k = k+1;   toc;    
% end
% vfwiai = reshape(real(1./sqrt(fwiai)),n);

% k = 1; fwiqmk = 1./(v0(:)).^2;
% for f = f0:df:nf
%     tic;
%     model.f = f;   
%     fprintf('FWI Inversion frequency: %3d \n', f);
%     % misfit
%     tmp = ['D_' num2str(k) ];
%     dobs = eval(tmp);   
%     fh = @(m)misfit_fwiqq(m,dobs,alpha,model,sigm);     % FWI with qq^* covariance
%     
%     % Simple BB iteration
%     [fwiqmk,hist] = BBiter(fh,fwiqmk,1e-20,5);   
% 	dlmwrite(['../input/hist_fwiq_' num2str(f) '.txt'],hist);	
%  
%     k = k+1;   toc;    
% end
% vfwiq = reshape(real(1./sqrt(fwiqmk)),n);

%% plot
figure;fig1 = imagesc(x,z,v0+dv);colormap(jet);colorbar; xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       c = colorbar;c.Label.String = 'Velocity [Km/s]';set(gca,'fontsize',18); axis image       
% figure;fig2 = imagesc(x,z,v0);colormap(jet);colorbar; xlabel('Distance [m]','fontsize',18);
%        ylabel('Depth [m]','fontsize',18); hold on
%        c = colorbar;c.Label.String = 'Velocity [Km/s]';set(gca,'fontsize',18); axis image     
figure;fig3 = imagesc(x,z,vfwi,[2,2.3]);colormap(jet);colorbar; xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       c = colorbar;c.Label.String = 'Velocity [Km/s]';set(gca,'fontsize',18); axis image
figure;fig4 = imagesc(x,z,vwri,[2,2.3]);colormap(jet);colorbar; xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       c = colorbar;c.Label.String = 'Velocity [Km/s]';set(gca,'fontsize',18); axis image
% figure;fig5 = imagesc(x,z,vfwiai,[2,2.3]);colormap(jet);colorbar; xlabel('Distance [m]','fontsize',18);
%        ylabel('Depth [m]','fontsize',18); hold on
%        c = colorbar;c.Label.String = 'Velocity [Km/s]';set(gca,'fontsize',18); axis image
% figure;fig6 = imagesc(x,z,vfwiq,[2,2.3]);colormap(jet);colorbar; xlabel('Distance [m]','fontsize',18);
%        ylabel('Depth [m]','fontsize',18); hold on
%        c = colorbar;c.Label.String = 'Velocity [Km/s]';set(gca,'fontsize',18); axis image
       
%% save plots
% print(1,'-depsc','-r300',['../Fig/sim-vt']);
% print(2,'-depsc','-r300',['../Fig/v0']);
% print(3,'-depsc','-r300',['../Fig/sim-vi_fwi-ai']);
