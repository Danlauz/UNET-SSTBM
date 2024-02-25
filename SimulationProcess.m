function [Z,FV,ObsData]=SimulationProcess(x0,model,c,nbl,nbsim,ObsData,seed)

rng(seed)

%% Compute the spectral density 
%Spectral density computation
[F1,s,rot,cx]=DensSpec1Ddl(x0,model,c);

%% Frequency Vector
%1-Sampled isotropic random frequencies
UFV=rand(nbl,nbsim);
rf=interp1(F1{1},s{1},UFV);

%2-Quasi-randomly oriented vector over the unit-half sphere
rl=VanCorput(nbl,0); % Van Corput sequence

%3- Dilatation and roation to the anisotropic grid
DRGrid=cx{1}(:,[1 2 3])*rot{1};

%4- Random phases
Uphase=rand(nbl,nbsim);

FV.rf=rf;
FV.rl=rl;
FV.DRGrid=DRGrid;
FV.Uphase=Uphase;
FV.c=c;
FV.UFV=UFV;

%% Simulation process
Z=STBM(x0,FV,ObsData);

function Z=STBM(x0,FV,ObsData)
% Generated Gaussian random field using the spectral turning band method
%input:
%x0  -- Grid (nbpoint x dim)

% FV (structure containing frequency vectors info)
%FV.rf --   Random frequency of the isotropic covariance  (nl x nbsim)
%FV.rl -- Quasi random line orientation over the unit-half sphere (nl x 3)
%FV.DRGrid -- Dilation and rotation factor of grid points (n x 3)
%FV.U -- random uniform values (nl x nbsim)
%FV.c -- Sill

% ObsData (structure containing hard data info)
%ObsData.HD --   values of hard data  (nbHD)
%ObsData.LocHD -- localization of hard data (nbHD)
%ObsData.ki -- inverse kriging matrix of observation (nbHD x nbHD)
%ObsData.k0 -- kriging matrix ob simulated point (nbHD x n)

%return:
%Z  -- Random Gaussian field with prescribed covariance

% Authors : Dany Lauzon

[nl, nbsim]=size(FV.rf);

%nbsim simuation
for j=1:nbsim

    %Initialization
    zsim=zeros(size(x0,1),1);

    %Computation of each cosine function
    for i=1:nl
        zsim=zsim+ sqrt(FV.c)*cos(FV.DRGrid*(FV.rl(i,:).*FV.rf(i,j)')'+ 2*pi*FV.Uphase(i,j));
    end

    %Simulation zsim
    zsim=zsim*sqrt(2/nl);


    % Return Z
    Z(:,j)=zsim;
end

% Post-conditionning
if ~isempty(ObsData)
    Z= postcond( [x0(ObsData.LocHD,:) ,ObsData.HD] , [x0(ObsData.LocHD,:) , Z(ObsData.LocHD,:)], [x0 ,Z] , nbsim , ObsData.ki, ObsData.k0 );
    Z=Z(:,4:end);
end