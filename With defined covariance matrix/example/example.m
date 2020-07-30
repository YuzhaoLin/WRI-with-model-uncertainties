%% setup


% read model
dx = 10;  n = [81, 121];
v0 = 2*ones(n);              % initial model
dv = zeros(n);
dv(31:40,56:65) = 0.1;       % abnormal
dv(66:69,:)     = 0.14;     

m  = 1./(v0(:) + dv(:)).^2;
mk = 1./(v0(:)).^2;

% set frequency, do not set larger
%  than min(1e3*v(:))/(7.5*dx) 
%  or smaller than 0.5
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
model.zr = zr;   model.xr  = xr;    model.nw  = nf;
model.zs = zs;   model.xsf = xs;   

% generate matrices
P  = getP(model.h,model.n,model.zr,model.xr);

%% Noise in source and data
% data covariance
mum    = zeros(1,length(xr)); 
sigmm  = 2e-5*eye(length(xr),length(xr));
Mnoise = mvnrnd(mum, sigmm, length(xs))';

% source covariance
mus    = zeros(1,n(1)*n(2)); 
sigmp  = 2e-15*eye(n(1)*n(2),n(1)*n(2));
Pnoise = mvnrnd(mus, sigmp, length(xs))';

%% forward
for f = f0:df:nf
    model.f = f;
    fprintf('Forward frequency: %3d \n', f);
    d      = F(m,model,nf,Pnoise); 
    d      = d + Mnoise;
    tmp    = ['D_' num2str(k) ];
    eval([tmp,'=d;']);  k = k+1;
end

%% inversion
vel_tmp=zeros(n(1)*n(2),k-1);  ki = 1; 
for f = f0:df:nf
    tic;
    model.f = f;  
    fprintf('Inversion frequency: %3d \n', f);
    % misfit
    tmp = ['D_' num2str(ki) ];
    dobs = eval(tmp);
    fh = @(m)misfit(m,dobs,alpha,model);
   
    % Simple BB iteration
    [mk,hist] = BBiter(fh,mk,1e-30,5);     
    dlmwrite(['../input/hist_' num2str(f) '.txt'],hist);
    
    % modify update direction
    for icor = 1:n(1)*n(2) 
        if( (mk(icor))>max(m) )
           mk(icor) = max(m);
        end
        if( (mk(icor))<min(m) )
           mk(icor) = min(m);
        end
    end 
    vel_tmp(:,ki) = mk;
    ki = ki+1;   toc;    
end

vk = reshape(real(1./sqrt(mk)),n);

%% plot
figure;fig1 = imagesc(x,z,v0+dv);colormap(jet);colorbar; xlabel('Distance/m','fontsize',18);
       ylabel('Depth/m','fontsize',18); hold on
       c = colorbar;c.Label.String = 'Velocity(Km/s)';set(gca,'fontsize',18); axis image
       saveas(fig1,'../Fig/vt.fig');
figure;fig2 = imagesc(x,z,v0);colormap(jet);colorbar; xlabel('Distance/m','fontsize',18);
       ylabel('Depth/m','fontsize',18);  hold on
       c = colorbar;c.Label.String = 'Velocity(Km/s)';set(gca,'fontsize',18); axis image
       saveas(fig2,'../Fig/v0.fig');
figure;fig3 = imagesc(x,z,vk);colormap(jet);colorbar; xlabel('Distance/m','fontsize',18);
       ylabel('Depth/m','fontsize',18); hold on
       c = colorbar;c.Label.String = 'Velocity(Km/s)';set(gca,'fontsize',18); axis image
       saveas(fig3,'../Fig/v_fwi_un.fig');
cor = 60;
figure;fig4 = plot(v0(:,cor)+dv(:,cor),z,'LineWidth',2); hold on
        plot(v0(:,cor),z,'r','LineWidth',2); hold on
        plot(vk(:,cor),z,'g','LineWidth',2); set(gca,'YDir','reverse');
        xlabel('Velocity(m/s)','fontsize',18);ylabel('Depth(m)','fontsize',18); 
        saveas(fig4,'../Fig/vel_60.fig');
figure;for f = f0:df:nf
           tmp = ['../input/hist_' num2str(f) '.txt' ];
           hist = load(tmp);
           if(f == f0)
               ymax = max(hist(:,3));
           end
           fig5 = subplot(1,nf-f0+1,f-f0+1);
           plot(hist(:,3),'LineWidth',2);axis([0,length(hist(:,3)),0,ymax]);
           tmp = ['Iteration at ' num2str(f) ' Hz'];
           xlabel(tmp,'fontsize',15); ylabel('Misfit function','fontsize',15); 
        end  
        saveas(fig5,'../Fig/misfit.fig');

%% save plots
print(1,'-depsc','-r300',['../Fig/vt']);
print(2,'-depsc','-r300',['../Fig/v0']);
print(3,'-depsc','-r300',['../Fig/v_exfwi_est']);
print(4,'-depsc','-r300',['../Fig/vel_200_l']);
print(5,'-depsc','-r300',['../Fig/misfitl']);
