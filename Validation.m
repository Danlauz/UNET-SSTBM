%% Simulation of nbsim conditional simulation
nx=128; ny=128; nz=1;
x0=grille3(1,nx,1,1,ny,1,1,nz,1);
% Covariance model
model=[4 30 30 30 0 0 45]; c=[0.25];

% Numbers of simulated lines and number of simulation
nbl=4000;
nbsim=100;
seed=1259;
% Hard Data values
load('Reference_Data.mat')
nbHD=0;
LocHD = randperm(size(x0,1),nbHD);
HD=Zref(LocHD,1);

% Covariance matrix for the post-conditioning
k=covardm(x0(LocHD,1:3),x0(LocHD,1:3),model,c);
ki=inv(k);
k0=covardm(x0(LocHD,1:3),x0,model,c);

ObsData.LocHD=LocHD;
ObsData.HD=HD;
ObsData.ki=ki;
ObsData.k0=k0;

% Simulation of Gaussian random fields. 
[Zsim,FV,ObsData]=SimulationProcess(x0,model,c,nbl,nbsim,ObsData,seed);

BC=rand(2,nbsim)/5-0.1;
BC=[BC(1,:);BC(2,:)+1];

cR=rand(1,nbsim)*2-9;
cR=10.^cR;

Q=rand(1,nbsim)*10+10;
Q=-Q*litre/minute;

b=3;
n=0.3; %Aquifer porosity
cW=4.4*10^-10; % Water compressibility
rhow=1000; % Water density
g=9.81; %gravity
Ss=rhow*g*(cR+n*cW);
S=Ss*b;

%% Compute drawdown using flow simulator
for i=1:nbsim
    [R_Flow(:,:,i)]=FlowSimulation(nx,ny,Zsim(1:nx*ny,i),BC(:,i),cR(i),Q(i));
end

%% Flow Parameters for figures

%Piezometer location
x0Piezo=[];
x0Piezo(:,1)=nx/2 + [8 ; -16 ; -8  ;  16 ; 4 ; -4 ;...
                0  ;   0 ;  0  ;  0 ; 0 ; 0  ;...
                24 ;  24 ; -24 ; -24  ;...
                0  ;   0 ;  24 ; -24  ]*2;

x0Piezo(:,2)=ny/2 + [0 ;   0 ;   0 ;  0;  0 ; 0  ;...
                8  ; -16 ;  -8 ;  16; 4 ; -4 ;...
                24 ; -24 ; -24 ;  24 ;...
                24 ; -24 ;   0 ;   0  ]*2;
%Piezometer Index
LocData=nx*(x0Piezo(:,2)-1)+x0Piezo(:,1);

% Pumpimg and delta time at iteration
times=[0 1 5 15 30 60 120 180 720]*60;
pump=-0.003*meter^3/second*ones(1,9);pump(1)=0;
% Initialization of constant Data
ConstantData{1}=nx;
ConstantData{2}=ny;
ConstantData{3}=Rref(:,2:end)+randn(size(Rref(:,2:end)))*0.05;
ConstantData{4}=LocData; clear x0Piezo LocData


%% Loading Trained Network and Compute drawdown using U-Net
load("Image500T100V\\TrainedNetwork_500epoch.mat")
TrainingTime/60/60

R_AI=zeros(nx,ny,size(Rref,2)-1,nbsim);
for i=1:nbsim
    tic
    im1=reshape(Zsim(1:nx*ny,i),[nx ny]);
    im2=reshape(repmat(1:-1/(nx-1):0,1,ny),[nx ny])*(BC(2,i)-BC(1,i))+BC(1,i);
    for j=1:size(Rref,2)-1
        if j==1
            im3=reshape(repmat(1:-1/(nx-1):0,1,ny),[nx ny])*(BC(2,i)-BC(1,i))+BC(1,i);
            im5=zeros(nx,ny);
        else
            im3=R_AI(:,:,j-1,i);
            im5=zeros(nx,ny); im5(nx/2,ny/2)=Q(i)*minute/litre;
        end
        im4=zeros(nx,ny); im4(nx/2,ny/2)=times(j);
        im6=zeros(nx,ny); im6(nx/2,ny/2)=log10(Ss(i));
        R_AI(:,:,j,i)=predict(trainedNetwork,cat(3,im1,im2,im3,im4,im5,im6));
        %if j==1
        %    im3=reshape(repmat(1:-1/(nx-1):0,1,ny),[nx ny])*(BC(2,i)-BC(1,i))+BC(1,i);
        %else
        %    im3=reshape(R_Flow(:,j,i),[nx ny]);
        %end
        %R_AI_Flow(:,:,j,i)=predict(trainedNetwork,cat(3,im1,im2,im3,im4,im5,im6));
        % figure(j); imagesc(R_AI(:,:,j,i)); colorbar(); caxis([-1 1]); colormap('jet')
    end
    t(i)=toc;
end
mean(t)

%%
% Statistics
for j=1:nbsim
    xxx= reshape(R_AI(:,:,:,j),[nx*ny size(R_AI,3)]);
    for i=1:9
        RMSE(j,i)=sqrt(mean((R_Flow(:,i+1,j)-xxx(:,i)).^2,'all'));
        MAE(j,i)=mean(abs( R_Flow(:,i+1,j)-xxx(:,i) ),'all');

        % Sum of squared residuals
        SSR = sum((R_Flow(:,i+1,j)-xxx(:,i)).^2,'all');
        % Total sum of squares
        TSS = sum(((xxx(:,i) - mean(xxx(:,i),'all')).^2),'all');
        % R squared
        Rsquared(j,i) = 1 - SSR/TSS;
    end
end
mean(MAE,'all')
mean(MAE)
std(MAE)
mean(RMSE,'all')
mean(RMSE)
std(RMSE)
mean(Rsquared,'all')
mean(Rsquared)
std(Rsquared)
[mean(MAE,'all'), mean(RMSE,'all'), mean(Rsquared,'all'), TrainingTime/60/60]
%% Figure Validation Sim4
sim=78; 
for i=1:9
    figure(i)
    xxx= reshape(R_AI(:,:,:,sim),[nx*ny size(R_AI,3)]);
    plot(reshape(R_Flow(:,i+1,sim),[],1),reshape(xxx(:,i),[],1),'o','MarkerSize',6,...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
    hold on
    plot(-3.:0.01:1.1,(-3.:0.01:1.1),'-r','LineWidth',4)
    xlim([-3. 1.1])
    ylim([-3. 1.1])
    xlabel('Flow simulator')
    ylabel('U-Net Prediction')
    title('Scatterplot')
    set(gca,'LineWidth',1,'FontSize',14)
    %name=['C:\Users\dany-\OneDrive\Bureau\FigureArticle7\FigureValidation\Sim4\' 'Sim4_t' num2str(i) '.eps'];
    %saveas(gcf,name,'epsc')
end
figure(12)
xxx= reshape(R_AI(:,:,:,sim),[nx*ny size(R_AI,3)]);
plot(reshape(R_Flow(:,2:end,sim),[],1),reshape(xxx,[],1),'o','MarkerSize',6,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','black')
hold on
plot(-3:0.01:1,(-3:0.01:1),'-r','LineWidth',2)
hold on
plot(-3:0.01:1,(-3:0.01:1)-0.1,'--r','LineWidth',1)
hold on
plot(-3:0.01:1,(-3:0.01:1)+0.1,'--r','LineWidth',1)
xlim([-3 1])
ylim([-3 1])
xlabel('Flow simulator')
ylabel('U-Net Prediction')
set(gca,'LineWidth',1,'FontSize',14)

for im=1:9
    xxx= reshape(R_AI(:,:,:,sim),[nx*ny size(R_AI,3)]);
    LocData=ConstantData{4};
    %Source Localisation and Index
    x0Q=[nx/2 ny/2];
    LocInj=nx*(x0Q(:,2)-1)+x0Q(:,1);

    if im<=2
        contligne=[0,0,1,1,2,2,3,3,4,4]/5;
    elseif im>=3 && im<4
        contligne=[-2,-2,0,0,1,1,2,2,4,4]/5;
    elseif im>=4 && im<6
        contligne=[-4,-4,-2,-2,-1,-1,0,0,2,2,4,4]/5;
    elseif im>=6
        contligne=[-4,-4,-2,-2,-1,-1,0,0,2,2,4,4]/5;
    end

    figure(100*im+1)
    x = linspace(1,nx,ny);
    y = linspace(1,nx,ny);
    [X,Y] = meshgrid(x,y);
    Z = reshape(R_Flow(:,im+1,sim),[nx ny])';

    imagesc(Z);
    hold on
    plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
    hold on
    plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
    hold on
    [c1,h]=contour(X,Y,Z,contligne,'-k','ShowText','on','LineWidth',2);
    h.LevelList=round(h.LevelList,2);
    clabel(c1,h,'FontSize',10,'FontWeight','bold')
    colormap(jet)
    colorbar()
    clim([-2 1])
    xticklabels = 0:10:100;
    xticks = linspace(1, nx, numel(xticklabels));
    yticklabels = 0:10:100;
    yticks = linspace(1, ny, numel(yticklabels));
    set(gca,'YDir','normal', 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels,'LineWidth',1,'FontSize',14)
    axis off
    title('Flow Simulator')
    %name=['C:\Users\dany-\OneDrive\Bureau\FigureArticle7\FigureValidation\Sim4\' 'Sim4_t' num2str(im) '_FS.eps'];
    %saveas(gcf,name,'epsc')

    figure(100*im+2)
    x = linspace(1,nx,ny);
    y = linspace(1,nx,ny);
    [X,Y] = meshgrid(x,y);
    Z = reshape(xxx(:,im),[nx ny])';

    imagesc(Z);
    hold on
    plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
    hold on
    plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
    hold on
    [c1,h]=contour(X,Y,Z,contligne,'-k','ShowText','on','LineWidth',2);
    h.LevelList=round(h.LevelList,2);
    clabel(c1,h,'FontSize',10,'FontWeight','bold')
    colormap(jet)
    colorbar()
    clim([-2 1])
    xticklabels = 0:10:100;
    xticks = linspace(1, nx, numel(xticklabels));
    yticklabels = 0:10:100;
    yticks = linspace(1, ny, numel(yticklabels));
    title('U-Net prediction')
    set(gca,'YDir','normal', 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels,'LineWidth',1,'FontSize',14)
    axis off
    %name=['C:\Users\dany-\OneDrive\Bureau\FigureArticle7\FigureValidation\Sim4\' 'Sim4_t' num2str(im) '_Unet.eps'];
    %saveas(gcf,name,'epsc')

    figure(100*im+3)
    x = linspace(1,nx,ny);
    y = linspace(1,nx,ny);
    [X,Y] = meshgrid(x,y);
    Z = abs(reshape(xxx(:,im),[nx ny])'-reshape(R_Flow(:,im+1,sim),[nx ny])');

    imagesc(Z);
    hold on
    plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
    hold on
    plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
    hold on
    [c1,h]=contour(X,Y,Z,contligne,'-k','ShowText','on','LineWidth',2);
    h.LevelList=round(h.LevelList,2);
    clabel(c1,h,'FontSize',10,'FontWeight','bold')
    colormap(jet)
    colorbar()
    clim([0 0.2])
    xticklabels = 0:10:100;
    xticks = linspace(1, nx, numel(xticklabels));
    yticklabels = 0:10:100;
    yticks = linspace(1, ny, numel(yticklabels));
    set(gca,'YDir','normal', 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels,'LineWidth',1,'FontSize',14)
    axis off
end

figure(1000)
imagesc(reshape(Zsim(1:nx*ny,sim)+log10(3*darcy)+7,[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap('jet')
colorbar()
clim([-6 -3])
xticklabels = 0:10:100;
xticks = linspace(1, nx, numel(xticklabels));
yticklabels = 0:10:100;
yticks = linspace(1, ny, numel(yticklabels));
set(gca,'YDir','normal', 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels,'LineWidth',1,'FontSize',14)
axis on;
%name=['C:\Users\dany-\OneDrive\Bureau\FigureArticle7\FigureValidation\Sim4\' 'Sim4_K.eps'];
%saveas(gcf,name)

%% Figure Validation Sim3
sim=55;
for im=1:9
    xxx= reshape(R_AI(:,:,:,sim),[nx*ny size(R_AI,3)]);
    LocData=ConstantData{4};
    %Source Localisation and Index
    x0Q=[nx/2 ny/2];
    LocInj=nx*(x0Q(:,2)-1)+x0Q(:,1);

    if im<=2
        contligne=[0,0,1,1,2,2,3,3,4,4]/5;
    elseif im>=3 && im<4
        contligne=[-2,-2,-1,-1,0,0,1,1,2,2,3,3,4,4]/5;
    elseif im>=4 && im<6
        contligne=[-2,-2,-1,-1,0,0,1,1,2,2,3,3,4,4]/5;
    elseif im>=6
        contligne=[-4,-4,-2,-2,-1,-1,0,0,2,2,4,4]/5;
    end

    figure(100*im+1)
    x = linspace(1,nx,ny);
    y = linspace(1,nx,ny);
    [X,Y] = meshgrid(x,y);
    Z = reshape(R_Flow(:,im+1,sim),[nx ny])';

    imagesc(Z);
    hold on
    plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
    hold on
    plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
    hold on
    [c1,h]=contour(X,Y,Z,contligne,'-k','ShowText','on','LineWidth',2);
    h.LevelList=round(h.LevelList,2);
    clabel(c1,h,'FontSize',10,'FontWeight','bold')
    colormap(jet)
    colorbar()
    clim([-1 1])
    xticklabels = 0:10:100;
    xticks = linspace(1, nx, numel(xticklabels));
    yticklabels = 0:10:100;
    yticks = linspace(1, ny, numel(yticklabels));
    set(gca,'YDir','normal', 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels,'LineWidth',1,'FontSize',14)
    axis off
    title('Flow Simulator')
    %name=['C:\Users\dany-\OneDrive\Bureau\FigureArticle7\FigureValidation\Sim4\' 'Sim3_t' num2str(im) '_FS.eps'];
    %saveas(gcf,name,'epsc')

    figure(100*im+2)
    x = linspace(1,nx,ny);
    y = linspace(1,nx,ny);
    [X,Y] = meshgrid(x,y);
    Z = reshape(xxx(:,im),[nx ny])';

    imagesc(Z);
    hold on
    plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
    hold on
    plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
    hold on
    [c1,h]=contour(X,Y,Z,contligne,'-k','ShowText','on','LineWidth',2);
    h.LevelList=round(h.LevelList,2);
    clabel(c1,h,'FontSize',10,'FontWeight','bold')
    colormap(jet)
    colorbar()
    clim([-1 1])
    xticklabels = 0:10:100;
    xticks = linspace(1, nx, numel(xticklabels));
    yticklabels = 0:10:100;
    yticks = linspace(1, ny, numel(yticklabels));
    title('U-Net prediction')
    set(gca,'YDir','normal', 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels,'LineWidth',1,'FontSize',14)
    axis off
    %name=['C:\Users\dany-\OneDrive\Bureau\FigureArticle7\FigureValidation\Sim4\' 'Sim3_t' num2str(im) '_Unet.eps'];
    %saveas(gcf,name,'epsc')

    figure(100*im+3)
    x = linspace(1,nx,ny);
    y = linspace(1,nx,ny);
    [X,Y] = meshgrid(x,y);
    Z = abs(reshape(xxx(:,im),[nx ny])'-reshape(R_Flow(:,im+1,sim),[nx ny])');

    imagesc(Z);
    hold on
    plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
    hold on
    plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
    hold on
    [c1,h]=contour(X,Y,Z,contligne,'-k','ShowText','on','LineWidth',2);
    h.LevelList=round(h.LevelList,2);
    clabel(c1,h,'FontSize',10,'FontWeight','bold')
    colormap(jet)
    colorbar()
    clim([0 0.2])
    xticklabels = 0:10:100;
    xticks = linspace(1, nx, numel(xticklabels));
    yticklabels = 0:10:100;
    yticks = linspace(1, ny, numel(yticklabels));
    set(gca,'YDir','normal', 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels,'LineWidth',1,'FontSize',14)
    axis off
end

figure(1000)
imagesc(reshape(Zsim(1:nx*ny,sim)+log10(3*darcy)+7,[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap('jet')
colorbar()
clim([-6 -3])
xticklabels = 0:10:100;
xticks = linspace(1, nx, numel(xticklabels));
yticklabels = 0:10:100;
yticks = linspace(1, ny, numel(yticklabels));
set(gca,'YDir','normal', 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels,'LineWidth',1,'FontSize',14)
axis on;
%name=['C:\Users\dany-\OneDrive\Bureau\FigureArticle7\FigureValidation\Sim4\' 'Sim3_K.eps'];
%saveas(gcf,name)

for i=1:9
    figure(i)
    xxx= reshape(R_AI(:,:,:,sim),[nx*ny size(R_AI,3)]);
    plot(reshape(R_Flow(:,i+1,sim),[],1),reshape(xxx(:,i),[],1),'o','MarkerSize',6,...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
    hold on
    plot(-2.:0.01:1.1,(-2.:0.01:1.1),'-r','LineWidth',4)
    xlim([-2. 1.1])
    ylim([-2. 1.1])
    xlabel('Flow simulator')
    ylabel('U-Net Prediction')
    title('Scatterplot')
    set(gca,'LineWidth',1,'FontSize',14)
    %name=['C:\Users\dany-\OneDrive\Bureau\FigureArticle7\FigureValidation\Sim4\' 'Sim3_t' num2str(i) '.eps'];
    %saveas(gcf,name,'epsc')
end
figure(12)
xxx= reshape(R_AI(:,:,:,sim),[nx*ny size(R_AI,3)]);
plot(reshape(R_Flow(:,2:end,sim),[],1),reshape(xxx,[],1),'o','MarkerSize',6,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','black')
hold on
plot(-2:0.01:1,(-2:0.01:1),'-r','LineWidth',2)
hold on
plot(-2:0.01:1,(-2:0.01:1)-0.1,'--r','LineWidth',1)
hold on
plot(-2:0.01:1,(-2:0.01:1)+0.1,'--r','LineWidth',1)
xlim([-2 1])
ylim([-2 1])
xlabel('Flow simulator')
ylabel('U-Net Prediction')
set(gca,'LineWidth',1,'FontSize',14)
