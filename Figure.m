%Piezometer location
x0Piezo=[];
x0Piezo(:,1)=nx/2 + [8 ; -16 ; -8  ;  16 ; 4 ; -4 ;...
                    0  ;   0 ;  0  ;  0 ; 0 ; 0  ;...
                    24 ;  24 ; -24 ; -24  ;...
                    0  ;   0 ;  24 ; -24 ]*2;

x0Piezo(:,2)=ny/2 + [0 ;   0 ;   0 ;  0;  0 ; 0  ;...
                    8  ; -16 ;  -8 ;  16; 4 ; -4 ;...
                    24 ; -24 ; -24 ;  24 ;...
                    24 ; -24 ;   0 ;   0 ]*2;
%Piezometer Index
LocData=nx*(x0Piezo(:,2)-1)+x0Piezo(:,1);

%Source Localisation and Index
x0Q=[ nx/2 ny/2];
LocInj=nx*(x0Q(:,2)-1)+x0Q(:,1);

%% Figure
contligne=[-1,-1,0,0,1,1,2,2,3,3,4,4]/5;
figure(100)
x = linspace(1,nx,ny);
y = linspace(1,nx,ny);
[X,Y] = meshgrid(x,y);
Z = reshape(Rref(:,1),[nx,ny])';

imagesc(reshape(Rref(:,1),[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
hold on
[c1,h]=contour(X,Y,Z,contligne,'-k','ShowText','on','LineWidth',2);
h.LevelList=round(h.LevelList,2);
clabel(c1,h,'FontSize',15,'FontWeight','bold')
colormap(jet)
colorbar()
caxis([-0.5 1])
xticklabels = 0:10:100;
xticks = linspace(1, nx, numel(xticklabels));
yticklabels = 0:10:100;
yticks = linspace(1, ny, numel(yticklabels));
set(gca,'YDir','normal', 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels,'LineWidth',1.5,'FontSize',14)
axis on;
%title('Reference pressure head field (m)')
%%

figure(1)
imagesc(reshape(Zref(1:nx*ny)+log10(3*darcy)+7,[nx,ny])');
hold on
plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
hold on
plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
colormap('jet')
colorbar()
clim([-6 -3])
xticklabels = 0:10:100;
xticks = linspace(1, nx, numel(xticklabels));
yticklabels = 0:10:100;
yticks = linspace(1, ny, numel(yticklabels));
set(gca,'YDir','normal', 'XTick', xticks, 'XTickLabel', xticklabels, 'YTick', yticks, 'YTickLabel', yticklabels,'LineWidth',1,'FontSize',14)
axis on;
name=['C:\Users\dany-\OneDrive\Bureau\FigureArticle7\FigureRef\LogK_Ref.eps'];
saveas(gcf,name,'epsc')
%%
for i=1:size(Rref,2)-1
    figure(1+i)
    if i<=2
        contligne=[0,0,1,1,2,2,3,3,4,4]/5;
    elseif i>=3 && i<4
        contligne=[-1,-1,0,0,0.5,0.5,1,1,2,2,4,4]/5;
    elseif i>=4 && i<6
        contligne=[-2,-2,-1,-1,0,0,2,2,4,4]/5;
    elseif i>=6
        contligne=[-2,-2,-1,-1,0,0,2,2,4,4]/5;
    end
    x = linspace(1,nx,ny);
    y = linspace(1,nx,ny);
    [X,Y] = meshgrid(x,y);
    Z = reshape(Rref(:,1+i),[nx,ny])';

    imagesc(reshape(Rref(:,1+i),[nx,ny])');
    hold on
    plot(x0(LocInj,1),x0(LocInj,2),'ok','markersize',8,'LineWidth',2)
    hold on
    plot(x0(LocData,1),x0(LocData,2),'xk','markersize',8,'LineWidth',2)
    hold on
    [c1,h]=contour(X,Y,Z,contligne,'-k','ShowText','on','LineWidth',2);
    h.LevelList=round(h.LevelList,2);
    clabel(c1,h,'FontSize',15,'FontWeight','bold')
    colormap(jet)
    colorbar()
    caxis([-1 1])
    set(gca,'YDir','normal', 'LineWidth',1,'FontSize',16)
    axis off;
    %title('Reference pressure head field (m)')
    name=['C:\Users\dany-\OneDrive\Bureau\FigureArticle7\FigureRef\t' num2str(i) '.eps'];
    saveas(gcf,name,'epsc')
end
%%
figure(25)
TStep= [0 1 5 15 30 60 120 180 720];
plot([0 TStep],Rref(LocData,:),'k')
