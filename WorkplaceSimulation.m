%% Grid Data
nx=128; ny=128; nz=1;
x0=grille3(1,nx,1,1,ny,1,1,nz,1);


%% Covariance model
model=[4 30 30 30 0 0 45]; c=[0.25];

%% Numbers of simulated lines and number of simulation
nbl=4000;
nbsim=125; %125;250;375;500;625
seed=15248568;
%% Hard Data values
load('Reference_Data.mat','Zref')
nbHD=0;
LocHD = randperm(size(x0,1),nbHD);
HD=Zref(LocHD,1);

%2-covariance matrix for the post-conditioning
k=covardm(x0(LocHD,1:3),x0(LocHD,1:3),model,c);
ki=inv(k);
k0=covardm(x0(LocHD,1:3),x0,model,c);

ObsData.LocHD=LocHD;
ObsData.HD=HD;
ObsData.ki=ki;
ObsData.k0=k0;

%% Simulation of Gaussian random fields. 
[Z,FV,ObsData]=SimulationProcess(x0,model,c,nbl,nbsim,ObsData,seed);

%% Figure

LocHD=ObsData.LocHD;

figure(1)
imagesc(reshape(Z(1:nx*ny,1),[nx ny])')
hold on 
plot(x0(LocHD,1),x0(LocHD,2),'ko')
colorbar()

BC=rand(2,nbsim)/5-0.1;
BC=[BC(1,:);BC(2,:)+1];

cR=rand(1,nbsim)*2-9;
cR=10.^cR;

Q=rand(1,nbsim)*10+10;
Q=-Q*litre/minute;

%%
for i=1:nbsim
    i
    R(:,:,i)=FlowSimulation(nx,ny,Z(:,i),BC(:,i),cR(i),Q(i));
end
%%
for i=1:nbsim
    Data(:,:,1,i)=reshape(Z(1:nx*ny,i),[nx ny]);
end
for i=1:nbsim
    for j=1:size(R,2)
        Data(:,:,j+1,i)=reshape(R(:,j,i),[nx ny]);
    end
end

%%
save('DataTrainValidation.mat','Data','BC','cR','Q')

