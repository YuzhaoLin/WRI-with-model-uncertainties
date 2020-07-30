function [data] = read(n1,n2,filename)

%%
% This function is to read .dat file
%    n1 = vertical cordinate
%    n2 = horizontal cordinate
%

%%

fp   =  fopen(filename,'rb');
data =  fread(fp,[n1,n2],'float32');
fclose(fp);

% figure;imagesc(10*x,10*z,a);colorbar;colormap(jet);
% xlabel('Distance/m','fontsize',10);
% ylabel('Depth/m','fontsize',10); title('FWI Result');
% set(gca,'fontsize',10);%set(gca,'XAxisLocation','top'); 

% figure; plot(vt(:,100),'LineWidth',2); hold on
%           plot(v0(:,100),'r','LineWidth',2); hold on
%           plot(vp(:,100),'g','LineWidth',2);

% figure; plot(v(:,100),(1:depth)*dz,'LineWidth',2);
% set(gca,'YDir','reverse');
% xlabel('Velocity(m/s)','fontsize',10);
% ylabel('Depth(m)','fontsize',10); title('Velocity Curve');
   

end

