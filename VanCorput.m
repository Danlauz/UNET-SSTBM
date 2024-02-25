function line=VanCorput(nbline,init)

if nargin==1
    init=0;
end

count=0;
for i= 1+init : nbline+init
    count=count+1;
    
    a2 = base(i,2);
    a3 = base(i,3);
    u(count,1) = sum(a2./(2.^[1 : size(a2,2)]));
    u(count,2) = sum(a3./(3.^[1 : size(a3,2)]));
end
line = [cos(2*pi*u(:,1)).*sqrt(1-u(:,2).^2),sin(2*pi*u(:,1)).*sqrt(1-u(:,2).^2),u(:,2)];


function y=base(n,b)
% trouver k
kmax=floor(log(n)/log(b));
n2=n;
for k=kmax:-1:0
   y(k+1)=floor(n2/b^k);
   n2=mod(n2,b^k);
end