% T/S plots for nutrients from the Kuroshio
% transects

% clear all
addpath /users/sclayton/Documents/MATLAB/matlab/seawater

%% load data
load ../data/ostreo.mat

O=log10(O(1:31,:));
O(O==-Inf)=0;
kuro=find(lon>130);

%% make a matrix of sigma_t
SS=33:0.1:34.7;
TT=13:0.05:28;
SS=repmat(SS,length(TT),1);
TT=repmat(TT,size(SS,2),1)';
sigt=sw_dens(SS,TT,0);
sigt=sigt-1.*1000;

%% plot the data
figure(1)

subplot(1,2,1)
[C]=contour(SS,TT,sigt,20,'k');
clabel(C);
hold on
scatter(S(kuro),theta(kuro),1000,O(:,1),'.');colorbar;caxis([0 4.7])
xlabel('Salinity')
ylabel('Potential Temperature')
title('OI (log_{10} cells ml^{-1})')
hold off

subplot(1,2,2)
[C]=contour(SS,TT,sigt,20,'k');
clabel(C);
hold on
scatter(S(kuro),theta(kuro),1000,O(:,2),'.');colorbar;caxis([0 4.7])
xlabel('Salinity')
ylabel('Potential Temperature')
title('OII (log_{10} cells ml^{-1})')
hold off

figure(2)
colormap default, cmap=colormap; colormap(cmap([15  20 24:2:58],:));
[C]=contour(SS,TT,sigt,20,'k');
clabel(C);
hold on
scatter(S(kuro),theta(kuro),4000,ostreo(kuro)./100,'.');
colorbar;caxis([0 1])
xlabel('Salinity','FontSize',16)
ylabel('Temperature','FontSize',16)
%title('%OII clade')
set(gca,'FontSize',16,'FontName','Helvetica')
hold off
