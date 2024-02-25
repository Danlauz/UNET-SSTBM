%% Grid Data
nx=128; ny=128; nz=1;
x0=grille3(1,nx,1,1,ny,1,1,nz,1);
seed=45816;
%% Covariance model
model=[4 30 30 30 0 0 0]; c=[0.25];

%% Numbers of simulated lines and number of simulation
nbl=3000;
nbsim=1;

%% Hard Data values
x0Q=[nx/2 ny/2]; LocHD=nx*(x0Q(:,2)-1)+x0Q(:,1);
HD=[-0.25];

%2-covariance matrix for the post-conditioning
k=covardm(x0(LocHD,1:3),x0(LocHD,1:3),model,c);
ki=inv(k);
k0=covardm(x0(LocHD,1:3),x0,model,c);

ObsData.LocHD=LocHD;
ObsData.HD=HD;
ObsData.ki=ki;
ObsData.k0=k0;

%% Simulation of Gaussian random fields. 
[Zref,FVref]=SimulationProcess(x0,model,c,nbl,nbsim,ObsData,seed);

%% Figure
figure(1)
imagesc(reshape(Zref(1:nx*ny,1),[nx ny])')
colorbar()

BC=[-0.03509; 1.0423];
cR= 8.1623e-09;
Q=-17*litre/minute;
%%
Rref=FlowSimulation(nx,ny,Zref,BC,cR,Q);

%%

save('Reference_Data.mat','Zref','Rref')

Figure