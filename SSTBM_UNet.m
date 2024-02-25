function [Z,Err,Iter]=SSTBM_UNet(x0,model,c,nbsim,nl,niter,seed,PC,ConstantData)
% Generated calibrated Gaussian random field using the sequential spectral turning band method
    %input:
    %x0  -- Grid (nbpoint x dim)
    %model -- covariance model
    %c -- sill (one for each model)
    %nbsim -- number of simulation
    %nl -- number of lines calibrated
    %MaxIter -- maximum number of iteration for the optimizer
    %type -- Objective function to evaluate
    %OFmin -- minimum objective function value to reach
    %seed -- for reproductability
    %ConstantData -- Data needs for ObjectiveFunction.m

    %return:
    %Z  -- Calibrated Random Gaussian field
    %Err  -- Objective function value
    %Delta  -- Mean pixel perturbation
    %Phi  -- vector of optimized phase
    %Iter  -- Number of iteration performed
    
    % Authors : Dany Lauzon
    
%1- Spectral density computation
[SD1.F1,SD1.s,SD1.rot,SD1.cx]=DensSpec1Ddl(x0,model,c); SD1.c=c;


%2- Post-conditioning by simple kriging : Initialization
k=covardm(x0(PC.LocHD,:),x0(PC.LocHD,:),model,c);
PC.ki=inv(k);
PC.k0=covardm(x0(PC.LocHD,:),x0,model,c);

%nbsim simuation
for j=1:nbsim
    j       
    %For reproductability
    rng(475*j+seed);

    %Initialization
    SD=SD1;SD.U=[]; SD.z1{1}=[];
    err=nan(niter,1);
    errNow=1000;
    
    zsim=zeros(size(SD.cx{1},1),1);

    
    %Optimization on phase U
    for i=1:niter
        %3-Sampled random frequencies
        p=rand(nl,1);
        ul1=interp1(SD.F1{1},SD.s{1},p); % interpolate ul from p and cum
        %4-Random vector
        z=VanCorput(nl,(i-1)*nl+100);    % Van Corput sequence
        %5- Frequency vector
        z1=z.*ul1;

        Y=randn(2,nl);
        %6- Calibration process on the phase U
        options=optimset('MaxIter',2,'TolX',10^-8,'Display','off');
        func = @(t) OptErr(t,Y,zsim,z1,SD,ConstantData,x0,PC);
        X=rand(1);X=[X-0.2 X+0.2];
        [t,errNew] = fminbnd(func, X(1),X(2),options);
        U=normcdf(Y(1,:)*cos(2*pi*t)+Y(2,:)*sin(2*pi*t));

        if errNew<errNow
            errNow=errNew;
            SD.U=[SD.U U];
            SD.z1{1}=[SD.z1{1}; z1];
            zsim=zsim+sqrt(2)*sqrt(SD.c)* sum(cos((SD.cx{1}(:,[1 2 3])*SD.rot{1}')*z1'+ repmat(2*pi*U,size(x0,1),1)),2) ;
        end
        err(i)=errNow;

    end
    %7-Field Computation
    nlFinal=length(SD.U);
    zsim=zsim*sqrt(1/nlFinal);

    %8-If required : post-conditioning by kriging
    if ~isempty(PC.LocHD)
        z= postcond( [x0(PC.LocHD,:),PC.HD] , [x0(PC.LocHD,:) , zsim(PC.LocHD)], [x0 ,zsim] , 1 , PC.ki ,PC.k0 );
        zsim=z(:,end);
    end

    %9- Return Z and other parameters
    Err(:,j)=err;
    Z(:,j)=zsim;
    Iter(j)=i;
end

function [error]=OptErr(t,Y,zsim,z1,SD,ConstantData,x0,PC)

U=normcdf(Y(1,:)*cos(2*pi*t)+Y(2,:)*sin(2*pi*t));
nl=length([SD.U U]);

zsim=zsim+sqrt(2)*sqrt(SD.c)* sum(cos((SD.cx{1}(:,[1 2 3])*SD.rot{1}')*z1'+ repmat(2*pi*U,size(x0,1),1) ),2) ;

zsim=zsim*sqrt(1/nl);

if ~isempty(PC.LocHD)
    z= postcond( [x0(PC.LocHD,:),PC.HD] , [x0(PC.LocHD,:) , zsim(PC.LocHD)], [x0 ,zsim] , 1 , PC.ki ,PC.k0 );
    zsim=z(:,end);
end

error=ObjFunc_AI(zsim,ConstantData);

function [error,R_predict]=ObjFunc_AI(Z_sim,ConstantData)

nx=ConstantData{1};
ny=ConstantData{2};
R_ref=ConstantData{3};
LocData=ConstantData{4};
trainedNetwork=ConstantData{5};
times=ConstantData{6};
pump=ConstantData{7};

R_predict=zeros(nx,ny,size(R_ref,2));

BC=normcdf([Z_sim(nx*ny+1); Z_sim(nx*ny+2)]/sqrt(0.25))*0.2-0.1;
BC=[BC(1,:); BC(2,:)+1];

cR=normcdf([Z_sim(nx*ny+3)]/sqrt(0.25))*2-9;
cR=10.^cR;
b=3; n=0.3; %Aquifer porosity
cW=4.4*10^-10; % Water compressibility
rhow=1000; % Water density
g=9.81; %gravity
Ss=rhow*g*(cR+n*cW);

Q=normcdf([Z_sim(nx*ny+4)]/sqrt(0.25))*10+-20;

im1=reshape(Z_sim(1:nx*ny),[nx ny]);
im2=reshape(repmat(1:-1/(nx-1):0,1,ny),[nx ny])*(BC(2)-BC(1))+BC(1);
for j=1:size(R_ref,2)
    if j==1
        im3=reshape(repmat(1:-1/(nx-1):0,1,ny),[nx ny])*(BC(2)-BC(1))+BC(1);
        im5=zeros(nx,ny);
    else
        im3=R_predict(:,:,j-1);
        im5=zeros(nx,ny); im5(nx/2,ny/2)=Q;
    end
    im4=zeros(nx,ny); im4(nx/2,ny/2)=times(j);
    im6=zeros(nx,ny); im6(nx/2,ny/2)=log10(Ss);

    R_predict(:,:,j)=predict(trainedNetwork,cat(3,im1,im2,im3,im4,im5,im6));
end

R_predict=reshape(R_predict,[nx*ny size(R_ref,2)]);

error= sqrt( mean( (R_ref(LocData,:)-R_predict(LocData,:)).^2 ,'all') ) ;
