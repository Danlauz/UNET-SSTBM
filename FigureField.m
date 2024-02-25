
%% Figure U-Net prediction 
%Source Localisation and Index
x0Q=[nx/2 ny/2];
LocInj=nx*(x0Q(:,2)-1)+x0Q(:,1);

figure(201)
imagesc(reshape(Zref-4.5286,[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap(jet)
colorbar()
xticklabels = 0:10:100;
xticks = linspace(1, 101, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:10:100;
yticks = linspace(1, 101, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca,'YDir','normal')
clim([-6 -3])
axis off;
set(gca,'FontSize',14)

sim=34;
figure(202)
imagesc(reshape(Z_AI(1:nx*ny,sim)-4.5286,[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap(jet)
colorbar()
xticklabels = 0:10:100;
xticks = linspace(1, 101, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:10:100;
yticks = linspace(1, 101, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca,'YDir','normal')
title('U-Net: Simulation #1')
clim([-6 -3])
axis off;
set(gca,'FontSize',14)

figure(203)
imagesc(reshape(ZPredict_Flow(1:nx*ny,sim)-4.5286,[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap(jet)
colorbar()
xticklabels = 0:10:100;
xticks = linspace(1, 101, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:10:100;
yticks = linspace(1, 101, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca,'YDir','normal')
title('MRST: Simulation #1')
clim([-6 -3])
axis off;
set(gca,'FontSize',14)

sim=87;
figure(204)
imagesc(reshape(Z_AI(1:nx*ny,sim)-4.5286,[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap(jet)
colorbar()
xticklabels = 0:10:100;
xticks = linspace(1, 101, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:10:100;
yticks = linspace(1, 101, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca,'YDir','normal')
title('U-Net: Simulation #2')
clim([-6 -3])
axis off;
set(gca,'FontSize',14)

figure(205)
imagesc(reshape(ZPredict_Flow(1:nx*ny,sim)-4.5286,[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap(jet)
colorbar()
xticklabels = 0:10:100;
xticks = linspace(1, 101, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:10:100;
yticks = linspace(1, 101, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca,'YDir','normal')
title('MRST: Simulation #2')
clim([-6 -3])
axis off;
set(gca,'FontSize',14)

figure(206)
imagesc(reshape(mean(Z_AI(1:nx*ny,:),2)-4.5286,[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap(jet)
colorbar()
xticklabels = 0:10:100;
xticks = linspace(1, 101, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:10:100;
yticks = linspace(1, 101, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca,'YDir','normal')
title('Mean logK using U-Net')
clim([-6 -3])
axis off;
set(gca,'FontSize',14)

figure(207)
imagesc(reshape(mean(ZPredict_Flow(1:nx*ny,:),2)-4.5286,[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap(jet)
colorbar()
xticklabels = 0:10:100;
xticks = linspace(1, 101, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:10:100;
yticks = linspace(1, 101, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca,'YDir','normal')
title('Mean logK using MRST')
clim([-6 -3])
axis off;
set(gca,'FontSize',14)


figure(306)
imagesc(reshape(var(Z_AI(1:nx*ny,:),1,2),[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap(jet)
%colorbar()
xticklabels = 0:10:100;
xticks = linspace(1, 101, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:10:100;
yticks = linspace(1, 101, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca,'YDir','normal')
title('Variance U-Net')
axis off;
clim([0 0.30])
colorbar();
set(gca,'FontSize',14)

figure(308)
imagesc(reshape(var(ZPredict_Flow(1:nx*ny,:),1,2),[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap(jet)
%colorbar()
xticklabels = 0:10:100;
xticks = linspace(1, 101, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:10:100;
yticks = linspace(1, 101, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca,'YDir','normal')
title('Variance MRST')
axis off;
clim([0 0.30])
colorbar();
set(gca,'FontSize',14)
%%
for i=1:nbsim
   corr2_Unet(i)=corr2(Zref,Z_AI(:,i));
   corr2_MRST(i)=corr2(Zref,ZPredict_Flow(:,i));
end
mean(corr2_Unet)
std(corr2_Unet)
mean(corr2_MRST)
std(corr2_MRST)