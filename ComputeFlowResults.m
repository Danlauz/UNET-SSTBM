%% Calibration process using flow simulator
load('ZPredict_Flow_400T80V_200itt.mat')

load("Image400T80V\\TrainedNetwork_200epoch.mat")
load("ZPredict400T80V_200epoch.mat")
%% Parameters
load('Reference_Data.mat')
nx=128; ny=128;nz=1;
x0=grille3(1,nx,1,1,ny,1,1,nz,1);
seed=85632;

%Piezometer location
x0Piezo=[];
x0Piezo(:,1)=nx/2 + [8 ; -16 ; -8  ;  16 ; 4 ; -4 ;...
                0  ;   0 ;  0  ;  0 ; 0 ; 0  ;...
                24 ;  24 ; -24 ; -24  ;...
                0  ;   0 ;  24 ; -24  ]*2;

x0Piezo(:,2)=ny/2 + [0 ;   0 ;   0 ;  0;  0 ; 0  ;...
                8 ; -16 ;  -8 ;  16; 4 ; -4 ;...
                24 ; -24 ; -24 ;  24 ;...
                24 ; -24 ;   0 ;   0  ]*2;
%Piezometer Index
LocData=nx*(x0Piezo(:,2)-1)+x0Piezo(:,1);

% Covariance model
model=[4 30 30 30 0 0 45]; c=0.25;

% Number of simulation
nbsim=100;

RrefData=Rref;
RrefData=Rref(LocData,1:end)+error;

% Pumpimg and delta time at iteration
times=[0 1 5 15 30 60 120 180 720]*60;
pump=-0.003*meter^3/second*ones(1,9);pump(1)=0;
% Initialization of constant Data
ConstantData{1}=nx;
ConstantData{2}=ny;
ConstantData{3}=RrefData(:,2:end);
ConstantData{4}=LocData; clear x0Piezo LocData
ConstantData{5}=trainedNetwork; 
ConstantData{6}=times; clear times
ConstantData{7}=pump; clear pump

%% Compute drawdown using U-Net and calibrated Gaussian field
BC_AI=normcdf([Z_AI(nx*ny+1,:); Z_AI(nx*ny+2,:)]/sqrt(0.25))/5-0.1;
BC_AI=[BC_AI(1,:);BC_AI(2,:)+1];

Cr_AI=normcdf([Z_AI(nx*ny+3,:)]/sqrt(0.25))*2-9;
Cr_AI=10.^Cr_AI;

b=3; n=0.3; %Aquifer porosity
cW=4.4*10^-10; % Water compressibility
rhow=1000; % Water density
g=9.81; %gravity
Ss_AI=rhow*g*(Cr_AI+n*cW);

Q_AI=normcdf([Z_AI(nx*ny+4,:)]/sqrt(0.25))*10+-20;

BC_Predict_Flow=normcdf([ZPredict_Flow(nx*ny+1,:); ZPredict_Flow(nx*ny+2,:)]/sqrt(0.25))/5-0.1;
BC_Predict_Flow=[BC_Predict_Flow(1,:);BC_Predict_Flow(2,:)+1];

Cr_Predict_Flow=normcdf([ZPredict_Flow(nx*ny+3,:)]/sqrt(0.25))*2-9;
Cr_Predict_Flow=10.^Cr_Predict_Flow;

b=3; n=0.3; %Aquifer porosity
cW=4.4*10^-10; % Water compressibility
rhow=1000; % Water density
g=9.81; %gravity
Ss_Predict_Flow=rhow*g*(Cr_Predict_Flow+n*cW);

Q_Predict_Flow=normcdf([ZPredict_Flow(nx*ny+4,:)]/sqrt(0.25))*10+-20;


R_AI=zeros(nx,ny,size(Rref,2)-1,nbsim);
for i=1:nbsim
    tic
    im1=reshape(Z_AI(1:nx*ny,i),[nx ny]);
    im2=reshape(repmat(1:-1/(nx-1):0,1,ny),[nx ny])*(BC_AI(2,i)-BC_AI(1,i))+BC_AI(1,i);
    for j=1:size(Rref,2)-1
        if j==1
            im3=reshape(repmat(1:-1/(nx-1):0,1,ny),[nx ny])*(BC_AI(2,i)-BC_AI(1,i))+BC_AI(1,i);
            im5=zeros(nx,ny);
        else
            im3=R_AI(:,:,j-1,i);
            im5=zeros(nx,ny); im5(nx/2,ny/2)=Q_AI(i);
        end
        im4=zeros(nx,ny); im4(nx/2,ny/2)=ConstantData{6}(j);
        im6=zeros(nx,ny); im6(nx/2,ny/2)=log10(Ss_AI(i));
        R_AI(:,:,j,i)=predict(trainedNetwork,cat(3,im1,im2,im3,im4,im5,im6));
    end
    t=toc;
end
%% Compute drawdown using flow simulator and calibrated Gaussian field
for i=1:nbsim
   [R_AI_Flow(:,:,i)]=FlowSimulation(nx,ny,Z_AI(1:nx*ny,i),BC_AI(:,i),Cr_AI(i),Q_AI(i)*litre/minute);
   [R_Flow(:,:,i)]=FlowSimulation(nx,ny,ZPredict_Flow(1:nx*ny,i),BC_Predict_Flow(:,i),Cr_Predict_Flow(i),Q_Predict_Flow(i)*litre/minute);
end
%%
LocData=ConstantData{4};
sigma=0.05;
figure(1)
xxx= reshape(R_AI(:,:,:,:),[nx*ny size(R_AI,3) nbsim]);
plot(reshape(Rref(LocData,2:end),[],1),reshape(mean(xxx(LocData,:,:),3),[],1),'ob','MarkerSize',4,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','black')
hold on
x=reshape(Rref(LocData,2:end),[],1);
y=reshape(mean(xxx(LocData,:,:),3),[],1);
neg=y-reshape(quantile(xxx(LocData,:,:),0.05,3),[],1);
pos=reshape(quantile(xxx(LocData,:,:),0.95,3),[],1)-y;
errorbar(x,y,neg,pos,'+k')
hold on 
p1=plot(-0.85:0.01:1,(-0.85:0.01:1)+sigma*1.96,'--r','LineWidth',1);
hold on
plot(-0.85:0.01:1,(-0.85:0.01:1)-sigma*1.96,'--r','LineWidth',1)
hold on
plot(-0.85:0.01:1,(-0.85:0.01:1),'-r','LineWidth',1)
xlim([-0.85 1])
ylim([-0.85 1])
legend([p1],{'C.I. 95%'}, 'Location', 'northwest')
xlabel('Reference')
ylabel('Inversion using U-Net')
set(gca,'FontSize',14)
%%

x0Q=[ nx/2 ny/2];
LocInj=nx*(x0Q(:,2)-1)+x0Q(:,1);

figure(306)
imagesc(reshape(var(Z_AI(1:nx*ny,:),1,2),[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap(jet)
for i=1:length(error)
    txt = ['e=' num2str(round(error(i)*100,1)) ' cm'];
    text(x0(LocData(i),1),x0(LocData(i),2),txt)
end
%colorbar()
xticklabels = 0:10:100;
xticks = linspace(1, 101, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:10:100;
yticks = linspace(1, 101, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca,'YDir','normal')
set(gca,'LineWidth',1.5,'FontSize',14)
axis off;
clim([0 0.30])
colorbar();
set(gca,'FontSize',14)

%%
figure(8)
h1 = histogram(rand(100,1)/5-0.1);
hold on
h2 = histogram(BC_AI(1,:));
hold on
plot([-0.03509 -0.03509], [0 30],'-k',LineWidth=3)

figure(9)
h1 = histogram(rand(100,1)/5-0.1+1);
hold on
h2 = histogram(BC_AI(2,:));
hold on
plot([1.0423 1.0423], [0 30],'-k',LineWidth=3)

figure(10)
h1 = histogram(rand(100,1)*2-9);
hold on
h2 = histogram(log10(Cr_AI(1,:)));
hold on
plot([log10(8.1623e-09) log10(8.1623e-09)], [0 30],'-k',LineWidth=3)

figure(11)
h1 = histogram(rand(100,1)*10+10);
hold on
h2 = histogram(-Q_AI(1,:));
hold on
plot([17 17], [0 30],'-k',LineWidth=3)
%%
[mean(BC_AI(1,:)), std(BC_AI(1,:)),quantile(BC_AI(1,:),0.15),quantile(BC_AI(1,:),0.85),skewness(BC_AI(1,:)),kurtosis(BC_AI(1,:))]
[mean(BC_AI(2,:)), std(BC_AI(2,:)),quantile(BC_AI(2,:),0.15),quantile(BC_AI(2,:),0.85),skewness(BC_AI(2,:)),kurtosis(BC_AI(2,:))]
[mean(log10(Ss_AI)), std(log10(Ss_AI)),quantile(log10(Ss_AI),0.15),quantile(log10(Ss_AI),0.85),skewness(log10(Ss_AI)),kurtosis(log10(Ss_AI))]
[mean(Q_AI), std(Q_AI),quantile(Q_AI,0.15),quantile(Q_AI,0.85),skewness(Q_AI),kurtosis(Q_AI)]
%%
figure(12)
h1 = histogram(rand(100,1)/5-0.1);
hold on
h2 = histogram(BC_Predict_Flow(1,:));
hold on
plot([-0.03509 -0.03509], [0 30],'-k',LineWidth=3)

figure(13)
h1 = histogram(rand(100,1)/5-0.1+1);
hold on
h2 = histogram(BC_Predict_Flow(2,:));
hold on
plot([1.0423 1.0423], [0 30],'-k',LineWidth=3)

figure(14)
h1 = histogram(rand(100,1)*2-9);
hold on
h2 = histogram(log10(Cr_Predict_Flow(1,:)));
hold on
plot([log10(8.1623e-09) log10(8.1623e-09)], [0 30],'-k',LineWidth=3)

figure(15)
h1 = histogram(rand(100,1)*10+10);
hold on
h2 = histogram(-Q_Predict_Flow(1,:));
hold on
plot([17 17], [0 30],'-k',LineWidth=3)
%%
[mean(BC_Predict_Flow(1,:)), std(BC_Predict_Flow(1,:)),quantile(BC_Predict_Flow(1,:),0.15),quantile(BC_Predict_Flow(1,:),0.85),skewness(BC_Predict_Flow(1,:)),kurtosis(BC_Predict_Flow(1,:))]
[mean(BC_Predict_Flow(2,:)), std(BC_Predict_Flow(2,:)),quantile(BC_Predict_Flow(2,:),0.15),quantile(BC_Predict_Flow(2,:),0.85),skewness(BC_Predict_Flow(2,:)),kurtosis(BC_Predict_Flow(2,:))]
[mean(log10(Ss_Predict_Flow)), std(log10(Ss_Predict_Flow)),quantile(log10(Ss_Predict_Flow),0.15),quantile(log10(Ss_Predict_Flow),0.85),skewness(log10(Ss_Predict_Flow)),kurtosis(log10(Ss_Predict_Flow))]
[mean(Q_Predict_Flow), std(Q_Predict_Flow),quantile(Q_Predict_Flow,0.15),quantile(Q_Predict_Flow,0.85),skewness(Q_Predict_Flow),kurtosis(Q_Predict_Flow)]