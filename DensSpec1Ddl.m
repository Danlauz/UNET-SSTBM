function [F1,s,R,x0rot]=DensSpec1Ddl(x0,model,c)

% here we define the equations for the various 3D spectral density isotrope with C=1. Any new model
% can be added here. note: f1=(2*pi*s)^2*f3
k=[];
Gam={'1'                                                                                ; %1-nugget ok
     '1./(pi^2*(1+s.^2).^2)'                                                            ; %2-exponential
     'exp(-s.^2/4)/(8*pi^(3/2))'                                                        ; %3-gaussian  
     '3/4/pi*besselj(3/2,s/2).^2./s.^3'                                                 ; %4-spherical  
     'N/A'                                                                              ; %5-linear
     '((s<=0.004).*0.0014776+(s>0.004).*(210./(pi^2.*s.^10)).*(6.*s.*cos(s./2)+(s.^2-12).*sin(s./2)).^2)'; %6-modele cubique
     ''                                                                                 ; %7-spline plaque mince
     'N/A'                                                                              ; %8-modele gravimétrique (Cauchy avec b=0.5) 
     'besselk(0,s)/pi^(3/2)/gamma(1.5)/4'                                               ; %9-modele magnétique (Cauchy avec b=1.5)
     ''                                                                                 ; %10-effet de trou sinusoidal
     ''                                                                                 ; %11-effet de trou cosinusoidal
     'N/A'                                                                              ; %12-christakos
     '((s<=0.1).*0.0012062+(s>0.1).*(60*(24.*cos(s)-s.^2.*cos(s)+9.*s.*sin(s)+4.*s.^2-24)/pi^2./s.^8))'; %13-wendland_1
     '((s<=0.05).*0.0009951+(s>0.05).*(27720./(pi^2.*s.^14)).*((s.^3-60.*s).*cos(s./2)+(120-12.*s.^2).*sin(s./2)).^2)'; %14-Penta
     '(s.^2+1).^(-5/2)*gamma(5/2)/pi^(3/2)/gamma(1)'                                    ; %15-Matern nu=1
     '(s.^2+1).^(-3)*gamma(3)/pi^(3/2)/gamma(3/2)'                                      ; %16-Matern nu=3/2
     '((s<=0.04).*0.000818757+(s>0.04).*(6720*(8*s.*(s.^2-24)+9*(35-2*s.^2).*sin(s)+s.*(s.^2-123).*cos(s)))/pi^2./s.^11)'  ;%17-wendland_2
     '((s<=0.002).*0.0016886+(s>0.002).*(2*s-3*sin(s)+s.*cos(s))/pi^2./s.^5)' }         ; %17-wendland_0
 
  
freqMax=[ '10       '; %1-nugget
          '100      '; %2-exponential
          '7        '; %3-gaussian
          '100      '; %4-spherical
          '         '; %5-linear
          '40       '; %6-modele cubique
          '         '; %7-spline plaque mince
          '         '; %8-modele gravimétrique (Cauchy avec b=0.5) 
          '15       '; %9-modele magnétique (Cauchy avec b=1.5)
          '         '; %10-effet de trou sinusoidal
          '100      '; %11-effet de trou cosinusoidal
          '         '; %12-christakos
          '30       '; %13-wendland_1
          '40       '; %14-Penta
          '20       '; %15-Matern nu=1
          '20       '; %16-Matern nu=3/2
          '40       '; %17-wendland_2
          '50       '];%18-wendland_0         

% some constants are defined
      
[n,d]=size(x0); % d dimension de l'espace
[rp,p]=size(c);
r=rp/p;  % nombre de structures
cx=[x0(:,1:d)];
nm=size(model,2);

% ne pas permettre des portées de 0 en input pour éviter divisions par 0
if nm>2
   model(:,2:1+d)=max(model(:,2:1+d),100*eps);
else
   model(:,2)=max(model(:,2),100*eps);
end

% calculer les densités spectrales
 f1i=zeros(n*p,1);
 for i=1:r
     
     % calculation of matrix of reduced rotated distances H
     [x0roti,Ri]=trans(x0,model,i);
     x0rot{i}=x0roti; R{i}=Ri;
     
     % initialization os frequence vector
     si=10.^[-5:0.00005:log10(str2num(freqMax(model(i,1),:)))];
     % Calculation of F from f1 from f3
     f1i=str2func(['@(s)(4*pi^2*s.^2).*(' Gam{model(i,1)}(:,:) ')']);     
     F1i=zeros(length(si),1);
     for j=1:length(si)
         F1i(j)=integral(f1i,eps,si(j));
     end
     F1i=F1i/F1i(end);
     
     s{i}=[0;si'];
     F1{i}=[0;F1i];
     f1{i}=f1i;
 end

function [cx,rot]=trans(cx,model,im)
% function [cx,rot]=trans(cx,model,im);
%
% TRANS is called from COKRI2. It takes as input original coordinates and
%       return the rotated and reduced coordinates following specifications
%       described in model(im,:)
% Rotations are all performed anticlockwise with the observer located on the positive side of 
% the axis and looking toward the origin. In 3D, rotations are performed first along z,
% then along rotated y and then along twice rotated x.
% Author: D. Marcotte
% Version 2.1  97/aug/18

% some constants are defined

[n,d]=size(cx);
[m,p]=size(model);

% check for 1-D or isotropic model

if p-1>d

   % perform rotation counterclockwise

   if d==2
      ang=model(im,4); cang=cos(ang/180*pi); sang=sin(ang/180*pi);
      rot=[cang,-sang;sang,cang];
   else

      % rotation matrix in 3-D is computed around z, y and x in that order

      angz=model(im,7); cangz=cos(angz/180*pi); sangz=sin(angz/180*pi);
      angy=model(im,6); cangy=cos(angy/180*pi); sangy=sin(angy/180*pi);
      angx=model(im,5); cangx=cos(angx/180*pi); sangx=sin(angx/180*pi);
      rotz=[cangz,-sangz,0;sangz,cangz,0;0 0 1];
      roty=[cangy,0,sangy;0 1 0;-sangy,0,cangy];
      rotx=[1 0 0;0 cangx -sangx;0 sangx cangx];
      rot=rotz*roty*rotx;
   end

   % rotation is performed around z, y and x in that order, the other coordinates are left unchanged.

   dm=min(3,d);
   cx(:,1:dm)=cx(:,1:dm)*rot;
   t=[model(im,2:1+dm),ones(d-dm,1)];
   t=diag(t);
 else
   t=eye(d)*model(im,2);
end

% perform contractions or dilatations (reduced h)

  cx=cx/t;