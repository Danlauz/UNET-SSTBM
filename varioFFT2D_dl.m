function [gh,nh]=varioFFT2D_dl(c,z,icode2,categ,display)
%
% function [Xmesh,Ymesh,gh,nh]=varioFFT2D(c,z,icode,categ,display); 
%
% function to compute variograms, cross-variograms, covariograms,
% cross-covariograms and pseudo-cross-variograms in 1D or 2D for up to 3 variables.
% the data are on a (possibly incomplete) regular grid
% the program computes variograms in the frequency domain using 2D-FFT. 
%
% INPUT : 
%
% c        n by d     matrix of coordinates (regular grid)
% z        n by nvar  matrix of values for the variables. Each line is
%                     associated with the corresponding line vector of
%                     coordinates in the c matrix, and each column corresponds
%                     to a different variable.
%                     Missing values are indicated by NaN
%                     IMPORTANT: At a point, either all variables are observed or all are missing
% icode2               a code to indicate which function to compute
%                      =1 : variograms and cross-variograms;
%                      =2 : covariograms and cross-covariograms
%                      =3 : variograms and pseudo-cross-variograms
%                      =4 : a mean is computed for the whole field instead of according to the lags
%                      =5 : covariance non-centrée (bivariate probabilities, for categorical data)
%                      =6 : transiograms (for categorical data)
%                      =7 : non-ergodic transiograms
%                      =8 : asymmetry from Bardossy and Horning (2017)
%                      =9 : asymmetry from Bardossy and Horning no ecdf(2017)
%                      =10: Rank correlation (2017)
% categ    boolean    tells if the variables are categorical (set to 1) or not
%                     (set to 0, default)
% display  boolean    tells if a plot must be displayed (set to 1) or not
%                     (set to 0, default)%                     
% OUTPUT :
%
% gh       nvar by nvar cell array of nx by ny (direct- and cross-) maps for
%                     variables i and j depending on icode.
% nh       nx by ny matrix of number of pairs 
%                     available to compute the structural function.
%
% This program uses the functions FFT2, IFFT2, FFT2SHIFT and CONJ which are
% standard MATLAB functions.
%
% Author: Dimitri D'Or - Ephesia Consult - 2014/11/17
% Modified from variof2.m written by D. Marcotte, dmarcotte@mail.polymtl.ca 
%
% Reference :
%
% Marcotte D., 1996. Fast Variogram computation with FFT. Computers & Geoscience, 22, 10, 1175-1186.
%space1=memory;
icode=icode2(1);
if nargin<5
    display=0;
end
if nargin<4
    categ=0;
end

[n_nodes,nvar]=size(z);
if categ
    nvar=length(unique(z(~isnan(z))));
end

%%% Finding the parameters of the grid

minc=min(c);
maxc=max(c);
dc=zeros(1,size(c,2));
for i=1:size(c,2)
    dum1=diff(unique(c(:,i)));
    if isempty(dum1)
        dc(i)=1;
    else
        dc(i)=dum1(1);
    end
end

nc=((maxc-minc)./dc)+1;

%%% Reformatting the data

if categ
    idnan=isnan(z);
    zc=zeros(n_nodes,nvar);
    for i=1:nvar
        zc(z==i,i)=1;
        zc(idnan,i)=nan;       % to keep the nan in the computation
    end
    z=zc;
    clear zc;
end

Z=cell(nvar,1);
Zid=cell(nvar,1);

if length(nc)==1
    nc=[nc 1];
end

for i=1:nvar
    Z{i}=reshape(z(:,i),nc);
end
%%% Initialization

if icode==6
    prop=zeros(1,nvar);
    for i=1:nvar
        prop(i)=nansum(z(:,i))/(length(z(:,i))-sum(idnan));   
    end
end

gh=cell(nvar,nvar);
[n,p]=size(Z{1});          % dimensions of data matrix

% find the closest multiple of 8 to obtain a good compromise between
% speed (a power of 2) and memory required

nrows=2*n-1;
ncols=2*p-1;
nr2=ceil(nrows/8)*8;
nc2=ceil(ncols/8)*8;
nv=[nr2,nc2];

% form an indicator  matrix:  1's for all data values
%                             0's for missing values
% in data matrix, replace missing values by 0;

for i=1:nvar
    Zid{i}=~isnan(Z{i});                        % 1 for a data value; 0 for missing
    Z{i}(~Zid{i})=0;                            % missing replaced by 0
end

% Preparation

fx=cell(nvar,1);
fx2=cell(nvar,nvar);
fxid=cell(nvar,nvar);
if icode==4
    m=zeros(1,nvar);
    for i=1:nvar
        m(i)=sum(sum(Z{i}(Zid{i})))/sum(sum(Zid{i}));
        Z{i}(Zid{i})=Z{i}(Zid{i})-m(i);
    end
end

% Compute number of pairs

for i=1:nvar
    %for j=1:nvar
    j=i;
        if i==j
            fx{i}=fftn(Z{i},nv);                         % fourier transform of Z{i}
            fxid{i,i}=fftn(Zid{i},nv);
        else
            fxid{i,j}=fftn(Zid{i}.*Zid{j},nv);
        end
        fx2{i,j}=fftn(Z{i}.*Z{j},nv);                    % fourier transform of Z{i}*Z{j}
    %end
end

%if categ  % for future use to consider the case where not all variables are present at data point (will require more modifications elsewhere)
nh=round(real(ifftn(conj(fxid{1,1}).*fxid{1,1})));           % number of pairs
%else
%     for i=1:nvar
%         for j=1:nvar
%             nh{i,j}=round(real(ifftn(conj(fxid{i,j}).*fxid{i,j})));            % number of pairs
%         end
%     end
% end

% compute the different structural functions according to icode

switch icode
    case 1  % variograms and cross-variograms are computed        
        for i=1:nvar
            j=i;
            %for j=1:nvar
                t1=fftn(Z{i}.*Zid{j},nv);
                t2=fftn(Z{j}.*Zid{i},nv);
                t12=fftn(Z{i}.*Z{j},nv);
                gh{i,j}=real(ifftn(conj(fxid{i,j}).*t12+conj(t12).*fxid{i,j}-conj(t1).*t2-t1.*conj(t2)))./max(nh,1)/2;
            %end
        end
        
    case 2   % covariograms and cross-covariograms are computed
        for i=1:nvar
            j=i;
            %for j=1:nvar
            m_tail=real(ifftn(conj(fx{i}).*fxid{j,j}))./max(nh,1); % computes the tail means
            m_head=real(ifftn(conj(fxid{i,i}).*fx{j}))./max(nh,1); % computes the head means
            gh{i,j}=real((ifftn(conj(fx{i}).*fx{j}))./max(nh,1)-m_tail.*m_head);
            %end
        end
        
    case 3 % variograms and pseudo-cross-variograms are computed        
        for i=1:nvar
            for j=i:nvar
                gh{i,j}=real(ifftn(fxid{j,j}.*conj(fx2{i,i})+conj(fxid{i,i}).*fx2{j,j}-2*conj(fx{i}).*fx{j}))./max(nh,1)/2;
            end
        end
                
    case {4,5}
        for i=1:nvar
            j=i;
            %for j=1:nvar
            gh{i,j}=real((ifftn(conj(fx{i}).*fx{j}))./max(nh,1));
            %end
        end
                        
    case 6 % Transiograms are computed        
        for i=1:nvar
            for j=1:nvar
                gh{i,j}=(real((ifftn(conj(fx{i}).*fx{j}))./max(nh,1)))/prop(i);
            end
        end
        
    case 7  % non ergodic transiograms        
        fx_all=fx{1};
        for i=2:nvar
            fx_all=fx_all+fx{i};
        end
        
        for i=1:nvar
            propi=round(real(ifftn(conj(fx{i}).*fx_all)));
            for j=1:nvar
                gh{i,j}=real((ifftn(conj(fx{i}).*fx{j})))./max(propi,1);
            end
        end
        
    case 8  % Asymmetry Bardossy and Horning (2017)
        expo=icode2(2);
        Fz=cell(1,nvar);
        for i=1:nvar
            Fz{i}=phi_Z(z(:,i),nc);
            Fz{i}(~Zid{i})=0;
            sumPolNewton=0;
            for k=0:expo
                sumPolNewton= sumPolNewton+ ((-1).^k)*nchoosek(expo,k)*conj(fftn((Fz{i}.^(expo-k)).*Zid{i},nv)).*(fftn((Fz{i}.^k).*Zid{i},nv));
            end
            gh{i,i}=real(ifftn(sumPolNewton))./max(nh,1);
        end
        
    case 9  % Asymmetry Bardossy and Horning with Z{i} (2018)
        expo=3;%icode2(2);
        
        for i=1:nvar
            Z{i}(~Zid{i})=0;
            sumPolNewton=0;
            for k=0:expo
                sumPolNewton= sumPolNewton+ ((-1).^k)*nchoosek(expo,k)*conj(fftn((Z{i}.^(expo-k)).*Zid{i},nv)).*(fftn((Z{i}.^k).*Zid{i},nv));
            end
            gh{i,i}=real(ifftn(sumPolNewton))./max(nh,1);
        end
        
    case 10  % Rank correlation (2017)
        Fz=cell(nvar,1);
        for i=1:nvar
            Fz{i}=phi_Z(z(:,i),nc);
            Fz{i}(~Zid{i})=0;
            f1=fftn(Fz{i},nv);
            gh{i,i}=12*real(ifftn(conj(f1).*f1-0.5*conj(f1).*fxid{i,i}-0.5*conj(fxid{i,i}).*f1+0.25*conj(fxid{i,i}).*fxid{i,i}))./max(nh,1);
        end
        
    case 11  % High-order statistics, Cumulant tri-point (order 3)  (Article : Dimitrakopoulos (2009))
        
        clear nh
        for i=1:nvar
            % vector h1 and h2
            h1x=-1; h1y=0; %h1  %%%%north (-1,0)
            h2x=0; h2y=1; %h2   %%%% east (0,1)
            
            for factor=0:1:min(ceil(abs(n/h2x))-1,ceil(abs(p/h2y))-1)
                %finding fonction f(x).*f(x+h2)=g(x,x+h2) and pairs existence
                id=zeros(3*(n-1)+1,3*(p-1)+1);
                id(n:1:2*n-1,p:1:2*p-1)=Z{i};
                idh2=~(id(n:1:2*n-1,p:1:2*p-1).*id(n+factor*h2x:1:2*n-1+factor*h2x,p+factor*h2y:1:2*p-1+factor*h2y)==0);
                g=     id(n:1:2*n-1,p:1:2*p-1).*id(n+factor*h2x:1:2*n-1+factor*h2x,p+factor*h2y:1:2*p-1+factor*h2y);
                %some parameters
                t=nv/2+1;
                mid=((2*n-1)*(2*p-1)+1)/2;
                
                % finding number of tri-pairs (x,x+h1,x+h2) 
                if i==1
                    NH=round(real(ifftn(conj(fftn(idh2,nv)).*fxid{i,i})));
                    NH=fftshift(NH);                    
                    NH=NH(t(1)-n+1:t(1)+n-1,t(2)-p+1:t(2)+p-1);
                    num=    NH(mid:h1x+(2*p-1)*h1y: mid+(h1x+(2*p-1)*h1y)*min(ceil(abs(n/h1x))-1,ceil(abs(p/h1y))-1));
                    nh(:,factor+1)=num(1:end);
                end
                %  f(x)*f(x+h1)*f(x+h2)              
                GH=real((ifftn(conj(fftn(g,nv)).*fx{i})));                
                if isempty(GH)
                    continue;
                end
                ghtemp=fftshift(GH);
                ghtemp=ghtemp(t(1)-n+1:t(1)+n-1,t(2)-p+1:t(2)+p-1);               
                cum=ghtemp(mid:h1x+(2*p-1)*h1y: mid+(h1x+(2*p-1)*h1y)*min(ceil(abs(n/h1x))-1,ceil(abs(p/h1y))-1));               
                cum(num==0)=nan(1);                
                gh{i,i}(:,factor+1)=cum(1:end);
                
            end
            gh{i,i}(:,:)=gh{i,i}(:,:)./max(nh,1);
        end
        
        case 12  % High-order statistics, Cumulant tri-point (order 3)  (Article : Dimitrakopoulos (2009))
        clear nh
        for i=1:nvar
            %some parameters
            t=max(nv/2+1);
             fftZ1=fft(Z{i}(:,:),max(nv));
             fftID1=fft(Zid{i}(:,:),max(nv));
            
            % vector h1 and h2  %%% h1 direction : fix to north
            h2i=0; h2j=1; %h2   %%% h2 direction : all direction          
            for factor=0:1:min(ceil(abs(n/h2i))-1,ceil(abs(p/h2j))-1)
                %finding fonction f(x).*f(x+h2)=g(x,x+h2) and pairs existence
                id=zeros(3*(n-1)+1,3*(p-1)+1);
                id(n:1:2*n-1,p:1:2*p-1)=Z{i};   

                g=id(n:1:2*n-1,p:1:2*p-1).*id(n-factor*h2i:1:2*n-1-factor*h2i,p +factor*h2j:1:2*p-1+factor*h2j);
                idh2=~(g==0); 
                
                fftZ2=fft(g,max(nv));
                fftID2=fft(idh2,max(nv));
                
                NH=fftshift(ifft(sum(conj(fftID2).*fftID1,2)));
                GH=fftshift(ifft(sum(conj(fftZ2).*fftZ1,2)));
                nh(:,factor+1)=NH(t(1):-1:t(1)-n+1); 
                gh{i,i}(:,factor+1)=GH(t(1):-1:t(1)-n+1)./nh(:,factor+1);
            end   
        end
        
    case 13  % High-order statistics, Cumulant tri-point (order 3)  (Article : Dimitrakopoulos (2009))
        clear nh
        h2i=icode2(2); h2j=icode2(3); %h2   fix to one direction    
                
        Fzg=cell(1,nvar);
        Fzf=cell(1,nvar);
        for i=1:nvar                    
               
            
            %finding fonction f(x).*f(x+h2)=g(x,x+h2) and pairs existence
            id=zeros(3*(n-1)+1,3*(p-1)+1);
            id(n:1:2*n-1,p:1:2*p-1)=Z{i};
            g=id(n:1:2*n-1,p:1:2*p-1).*id(n+h2i:1:2*n-1+h2i,p+h2j:1:2*p-1+h2j);
            idh2=~(g==0);           
            
            %Transform to ecdf
            Fzg{i}=phi_Z(reshape(g,nc),nc);
            Fzf{i}=phi_Z(z(:,i),nc);
            Fzg{i}(~idh2)=0;
            Fzf{i}(~Zid{i})=0;
            % finding number of tri-pairs (x,x+h1,x+h2)
            if i==1
                nh=round(real(ifftn(conj(fftn(idh2,nv)).*fxid{i,i})));
            end
            %  f(x)*f(x+h1)*f(x+h2)
            gh{i,i}=real((ifftn(conj(fftn(Fzg{i},nv)).*fftn(Fzf{i},nv))))./max(nh,1);            
        end
        
        case 14  % Asymmetry Guthke and Bardossy (2016)
            for i=1:nvar,                
                Fz{i}=phi_Z(z(:,i),nc);
                Fz{i}(~Zid{i})=0;
                f3=fftn(Fz{i}.*Fz{i}.*Fz{i},nv);
                f2=fftn(Fz{i}.*Fz{i},nv);
                f1=fftn(Fz{i},nv);
                gh{i,i}=real(ifftn(conj(f3).*fxid{j,j} + 3.*conj(f2).*f1 - 3.*conj(f2).*fxid{j,j} + 3.*conj(f1).*f2 - 6.*conj(f1).*f1 + 3.*conj(f1).*fxid{j,j} + conj(fxid{i,i}).*f3 - 3.*conj(fxid{i,i}).*f2 + 3.*conj(fxid{i,i}).*f1-conj(fxid{i,i}).*fxid{j,j}))./max(nh,1);
                
            end;
        
        case 15  % Asymmetry Li (2010)
         for i=1:nvar
             Fz{i}=phi_Z(z(:,i),nc); 
             Fz{i}(~Zid{i})=0;
         end
        for i=1:nvar,
            for j=i:nvar,
                f2=fftn(Fz{i}.*Fz{i},nv);
                f1=fftn(Fz{i},nv);
                gh{i,i}=real(ifftn(conj(f2).*f1 - 0.5.*conj(f2).*fxid{j,j} - 2.*conj(f1).*f1 + conj(f1).*f2 + 0.75.*conj(f1).*fxid{j,j} -0.5.*conj(fxid{i,i}).*f2 + 0.75.*conj(fxid{i,i}).*f1 - 0.25.*conj(fxid{i,i}).*fxid{j,j} ))./max(nh,1);
            end;
        end;
        
end

% reduce matrices to required size
if icode ~= 11 && icode ~= 12
    t=nv/2+1;
    nh=fftshift(nh);
    nh=nh(t(1)-n+1:t(1)+n-1,t(2)-p+1:t(2)+p-1);
    for i=1:nvar
        for j=1:nvar
            if isempty(gh{i,j})
                continue;
            end
            ghtemp=fftshift(gh{i,j});
            gh{i,j}=ghtemp(t(1)-n+1:t(1)+n-1,t(2)-p+1:t(2)+p-1);
        end
    end
end

% % Display graph

if display && (length(size(gh{1,1}))<3) % if display and we are in 2D
    figure
    for i=1:nvar
        for j=1:nvar
            if isempty(gh{i,j})
                continue;
            end
            subplot(nvar,nvar,j+nvar*(i-1));
            %      pcolor(Xmesh,Ymesh,gh{i,j});
            imagesc(gh{i,j})
            title(['Var. ',num2str(i),' vs. Var. ',num2str(j)]);
            axis equal
            %    axis([0 length(Xmesh) 0 length(Ymesh)]);
            shading flat;
            if icode==6
                caxis([0 1]);
                %             else
                %                 caxis([0 gmax]);
            end
            colorbar
        end
    end
end
%space2=memory;
%espace=space2.MemUsedMATLAB-space1.MemUsedMATLAB