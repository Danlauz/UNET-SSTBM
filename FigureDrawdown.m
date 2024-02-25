%% Figure Article Sim 1
sigma=0.05;
LocData=ConstantData{4};
sim=6;
figure(1)
xxx= reshape(R_AI(:,:,:,sim),[nx*ny size(R_AI,3)]);
plot(reshape(Rref(LocData,2:end),[],1),reshape(xxx(LocData,:),[],1),'o','MarkerSize',4,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','black')
hold on
plot(-0.90:0.01:1,(-0.90:0.01:1)+sigma*1.96,'--r','LineWidth',1)
hold on
plot(-0.90:0.01:1,(-0.90:0.01:1)-sigma*1.96,'--r','LineWidth',1)
hold on
p1=plot(-0.90:0.01:1,(-0.90:0.01:1),'-r','LineWidth',1);
xlim([-0.90 1])
ylim([-0.90 1])
legend([p1],'I.C. 95%','Location', 'northwest')
set(gca,'FontSize',14)
xlabel('Reference')
ylabel('Inversion Using U-Net')

figure(2)
xxx= reshape(R_Flow(:,2:end,sim),[nx*ny size(R_AI,3)]);
plot(reshape(Rref(LocData,2:end),[],1),reshape(xxx(LocData,:),[],1),'o','MarkerSize',4,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','black')
hold on
plot(-0.90:0.01:1,(-0.90:0.01:1)+sigma*1.96,'--r','LineWidth',1)
hold on
plot(-0.90:0.01:1,(-0.90:0.01:1)-sigma*1.96,'--r','LineWidth',1)
hold on
p1=plot(-0.90:0.01:1,(-0.90:0.01:1),'-r','LineWidth',1);
xlim([-0.90 1])
ylim([-0.90 1])
legend([p1],'I.C. 95%','Location', 'northwest')
xlabel('Reference')
ylabel('Inversion Using MRST')
set(gca,'FontSize',14)


% Figure Article Sim 2
LocData=ConstantData{4};
sim=74;
figure(11)
xxx= reshape(R_AI(:,:,:,sim),[nx*ny size(R_AI,3)]);
plot(reshape(Rref(LocData,2:end),[],1),reshape(xxx(LocData,:),[],1),'o','MarkerSize',4,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','black')
hold on
plot(-0.90:0.01:1,(-0.90:0.01:1)+sigma*1.96,'--r','LineWidth',1)
hold on
plot(-0.90:0.01:1,(-0.90:0.01:1)-sigma*1.96,'--r','LineWidth',1)
hold on
p1=plot(-0.90:0.01:1,(-0.90:0.01:1),'-r','LineWidth',1);
xlim([-0.90 1])
ylim([-0.90 1])
legend([p1],'I.C. 95%','Location', 'northwest')
xlabel('Reference')
ylabel('Inversion Using U-Net')
set(gca,'FontSize',14)

figure(12)
xxx= reshape(R_Flow(:,2:end,sim),[nx*ny size(R_AI,3)]);
plot(reshape(Rref(LocData,2:end),[],1),reshape(xxx(LocData,:),[],1),'o','MarkerSize',4,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','black')
hold on
plot(-0.90:0.01:1,(-0.90:0.01:1)+sigma*1.96,'--r','LineWidth',1)
hold on
plot(-0.90:0.01:1,(-0.90:0.01:1)-sigma*1.96,'--r','LineWidth',1)
hold on
p1=plot(-0.90:0.01:1,(-0.90:0.01:1),'-r','LineWidth',1);
xlim([-0.90 1])
ylim([-0.90 1])
legend([p1],'I.C. 95%','Location', 'northwest')
xlabel('Reference')
ylabel('Inversion Using MRST')
set(gca,'FontSize',14)

figure(21)
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

figure(22)
xxx= reshape(R_Flow(:,2:end,:),[nx*ny size(R_AI,3) nbsim]);
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
ylabel('Inversion using MRST')
set(gca,'FontSize',14)