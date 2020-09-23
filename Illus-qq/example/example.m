%% setup


% read model
dx = 10;  n = [21, 21];
v0 = 2*ones(n);             
mk = 1./(v0(:)).^2;

% set frequency, do not set larger
%  than min(1e3*v(:))/(7.5*dx) 
%  or smaller than 0.5
nf = 12;    f0 = 3;   df = 1;

% receivers
xr = 20:1*dx:190;
zr = 2*dx*ones(1,length(xr));

% sources
xs = 100;   
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
model.zs = zs;   model.xsf = xs;   

%% Noise in source and data
l     = wgn(1,n(1)*n(2),10);

sigp   = 2e-5;
mus    = 0; 
sigmi  = sigp*eye(n(1)*n(2),n(1)*n(2));
Pnoise0 = sigmi * l';

Q       = getQ_for(model.h,model.n,model.zs,model.xsf,model.nf,5);
sigmq   = Q*Q';
Pnoiseq = sigmq * l';

Qe1 = zeros(n(1),n(2));
for iz = 1:n(1)
    for ix = 1:n(2)
        Qe1(iz,ix) = sqrt( (iz*dx-20)^2 + (ix*dx-100)^2 + 1);
    end
end
Qe1 = Qe1(:);
sigma1   = 1./(Qe1*Qe1');
Pnoisea1 = sigma1 * l';

Qe2 = zeros(n(1),n(2));
for iz = 1:n(1)
    for ix = 1:n(2)
        Qe2(iz,ix) =  sqrt( (iz*dx-20)^2 + (ix*dx-100)^2 );
    end
end
Qe2 = Qe2(:);
sigma2   = Qe2*Qe2';
Pnoisea2 = sigma2 * l';

%% plot
figure;imagesc(x,z,reshape(diag(sigmi),n));colormap(jet);colorbar; 
       xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       set(gca,'fontsize',18);  axis image   
figure;imagesc(x,z,reshape(diag(sigmq),n));colormap(jet);colorbar; 
       xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       set(gca,'fontsize',18);  axis image   
figure;imagesc(x,z,reshape(diag(sigma1),n));colormap(jet);colorbar; 
       xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       set(gca,'fontsize',18);  axis image   
figure;imagesc(x,z,reshape(diag(sigma2),n));colormap(jet);colorbar; 
       xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       set(gca,'fontsize',18);  axis image  
       
figure;imagesc(x,z,reshape(Pnoise0,n));colormap(jet);colorbar; 
       xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       set(gca,'fontsize',18);  axis image   
figure;imagesc(x,z,reshape(Pnoiseq,n));colormap(jet);colorbar; 
       xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       set(gca,'fontsize',18);  axis image   
figure;imagesc(x,z,reshape(Pnoisea1,n));colormap(jet);colorbar; 
       xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       set(gca,'fontsize',18);  axis image   
figure;imagesc(x,z,reshape(Pnoisea2,n));colormap(jet);colorbar; 
       xlabel('Distance [m]','fontsize',18);
       ylabel('Depth [m]','fontsize',18); hold on
       set(gca,'fontsize',18);  axis image   
        
%% save plots
print(1,'-depsc','-r300',['../Fig/sigp-i']);
print(2,'-depsc','-r300',['../Fig/sigp-q']);
print(3,'-depsc','-r300',['../Fig/sigp-a1']);
print(4,'-depsc','-r300',['../Fig/sigp-a2']);

print(5,'-depsc','-r300',['../Fig/noise-i']);
print(6,'-depsc','-r300',['../Fig/noise-q']);
print(7,'-depsc','-r300',['../Fig/noise-a1']);
print(8,'-depsc','-r300',['../Fig/noise-a2']);
