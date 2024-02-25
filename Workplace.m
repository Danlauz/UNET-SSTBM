%% Loading Data

% Grid Data
nx=128; ny=128;nz=1;
x0=grille3(1,nx,1,1,ny,1,1,nz,1);
x0=[x0 ;[-8525,-8945,1];[12142,98985,1];[-532684,7845963,1];[85496,-5415412,1]];
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

% Numbers of simulated lines and number of simulation
nl=100;
niter=400;
nbsim=100;
% Meeasurement errors
load('Reference_Data.mat')
rng(185425)
error=[-0.049;  0.045;  0.059;  0.015; -0.093;...
        0.096; -0.038; -0.027; -0.035; -0.033;...
       -0.043; -0.044; -0.011; -0.013;  0.025;...
        0.009; -0.021;  0.014;  0.019;  0.025];
RrefData=Rref;
RrefData(LocData,:)=RrefData(LocData,:)+error;

% Pumpimg and delta time at iteration
times=[0 1 5 15 30 60 120 180 720]*60;
pump=-0.003*meter^3/second*ones(1,9);pump(1)=0;

PC.LocHD=[];
PC.HD=[];


load('Image400T80V\\TrainedNetwork_200epoch.mat')

% Initialization of constant Data
ConstantData{1}=nx;
ConstantData{2}=ny;
ConstantData{3}=RrefData(:,2:end);
ConstantData{4}=LocData;
ConstantData{5}=trainedNetwork;
ConstantData{6}=times;
ConstantData{7}=pump; 
%Calibration process using U-Net prediction
tic
[Z_AI,Err_AI]=SSTBM_UNet(x0,model,c,nbsim,nl,niter,seed,PC,ConstantData);
t_AI=toc;
save('Zpredict400T80V_200epoch.mat','Z_AI','Err_AI','t_AI','error')

%%
tic
[ZPredict_Flow,Err_Flow]=SSTBM_Flow(x0,model,c,nbsim,nl,niter,seed,PC,ConstantData);
t_Flow=toc;
save('Zpredict_Flow_400T80V_200itt.mat','ZPredict_Flow','Err_Flow','t_Flow','error')

