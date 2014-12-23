% simple script to calculate concentration of ostreo across front using a
% simple reaction-diffusion model
% Sophie Clayton, October 2011, updated December 2014
% sclayton@mit.edu, sclayton@uw.edu

%clear all

n=24;

load ../data/ostreo % load the ostreo clade abundance data

kuro=find(lon>140 & z==0); % use only the data from the Kuroshio

% set up domain and set parameters
% SET UP GLOBAL PARAMETERS FOR THINGS THAT DON"T CHANGE BETWEEN MODEL RUNS

dy = 0.05;
y = 33.2:dy:34.4; % m
k =5e-08; % diffusity coeff in salinity space, S2/s
% k = 0.25*10^4; % diffusivity coeff, m2/s
% u1 =*10^-5.*ones(length(y)); % s-1
% u2 =u1;
%u2 = 1.1*10^-5.*exp(-((y-34).^2)./(2*1^2));
%u1 = 1*10^-5.*exp(-((y-33.8).^2)./(2*0.5^2));
%u1=0.9/(60*60*24).*ones(length(y));
u1=zeros(length(y));
u2=u1;
m1=0.1/(60*60*24);
%m1=0;
m2=m1;

par = [];

out = ostreo_model(par);

run('../figures/plot_ostreo_model.m')



%
%     nit=nitrate(ZZ==0 & stn<41);
%     
%    figure(2);
%    [AX,H1,H2] = plotyy(SALT,nit,y,u1,'plot'); set(AX,'FontSize',14,'XColor','k','YColor','k');
%    set(H1,'LineStyle','none','Marker','o','MarkerEdgeColor','k');
%    set(H2,'LineWidth',2,'Color','k')
%    % plot(y,u1,'k','LineWidth',2);axis([33.2 34.4 0 1.2*10^-5]);xlabel('salinity','FontSize',14);ylabel('\mu (day^{-1})','FontSize',14);set(gca,'FontSize',14);
% set(get(AX(1),'Ylabel'),'String','NO_3 + NO_2 (mmol m^{3})','FontSize',14,'Color','k')
% set(get(AX(2),'Ylabel'),'String','\mu (day^{-1})','FontSize',14,'Color','k')
% % end
