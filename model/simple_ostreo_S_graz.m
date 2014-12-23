% simple script to calculate concentration of ostreo across front
% Sophie Clayton, October 2011
% sclayton@mit.edu

clear all

load ../ostreo %/data1/sclayton/kuroshio/genomics/ostreo
%load /data1/sclayton/kuroshio/cruise_data/nitrate

kuro=find(lon>140 & z==0);

% set up domain and set parameters
dy = 0.05;
y = 33.2:dy:34.4; % m
k = 0.00000005; % S2/s
% k = 0.25*10^4; % m2/s
% u1 =*10^-5.*ones(length(y)); % s-1
% u2 =u1;
%u2 = 1.1*10^-5.*exp(-((y-34).^2)./(2*1^2));
u1 = 1*10^-5.*exp(-((y-33.8).^2)./(2*0.5^2));
%u1=0.9*10^-5.*ones(length(y));
%u1=zeros(length(y));
u2=u1;
k1=0.05;
m1=5*10^-6;
%m1=0;
m2=m1;
g=0.01;
mz=m1;



% time stepping
dt = 30; % quarter day in s
tend = 10;
tstep = 1:dt:60*60*24*tend;

OI=ones(length(y),length(tstep));
OII=ones(length(y),length(tstep));
Z=0.01.*ones(length(y),length(tstep));

% initial conditions and boundary conditions
OI(1,:)=100;
OII(end,:)=100;

% set gaussian nutrient profile
%k=10^4n0=5;
% n=0.5 + n0*exp(-((y-15000).^2)./(2*5000^2));
% umax = 1/(60*60*24);
% u1=(umax.*(n./(n+k1)));
% u2=u1;
n=ones(length(y),length(tstep));

for t = 2:length(tstep);
%    n0=5*(1+sin((2*pi*1/(60*60*24*3))*tstep(t)));
%    n(:,t)=0.5 + n0*exp(-((y-15000).^2)./(2*2000^2));
%    u1=(umax.*(n(:,t)./(n(:,t)+k1)));
%    u2=u1;
    for ystep = 2:length(y)-1;
        OI(ystep,t) = OI(ystep,t-1) + (k/(dy^2)*dt*(OI(ystep+1,t-1)-2*OI(ystep,t-1)+OI(ystep-1,t-1))+u1(ystep)*dt*OI(ystep,t-1))-m1*dt*OI(ystep,t-1)*Z(ystep,t-1);
        if OI(ystep,t)<0, OI(ystep,t)=0; end
        OII(ystep,t) = OII(ystep,t-1) + (k/(dy^2)*dt*(OII(ystep+1,t-1)-2*OII(ystep,t-1)+OII(ystep-1,t-1))+u2(ystep)*dt*OII(ystep,t-1))-m2*dt*OII(ystep,t-1)*Z(ystep,t-1);
        if OII(ystep,t)<0, OII(ystep,t)=0; end
        Z(ystep,t) = Z(ystep,t-1)*g*dt*(OI(ystep,t-1)+OII(ystep,t-1)) - mz*dt*Z(ystep,t-1) + k/(dy^2)*dt*(Z(ystep+1,t-1)-2*Z(ystep,t-1)+Z(ystep-1,t-1));
        if Z(ystep,t)<0, Z(ystep,t)=0; end
    end
end


figure(1);
% for t=length(tstep);
    subplot(3,1,1);plot(y,OI(:,end),'k',S(kuro),O(kuro,1),'ro','LineWidth',2);axis([33.2 34.4 0 18000]);title('OI','FontSize',14)
    ylabel('copies ml^{-1}','FontSize',14);set(gca,'FontSize',14)
    legend('model','data','North');%legend BOXOFF
%     hold on;plot(y,-(100/1.5)*y+2266.7,'r');hold off
    subplot(3,1,2);plot(y,OII(:,end),'k',S(kuro),O(kuro,2),'ro','LineWidth',2);axis([33.2 34.4 0 18000]);title('OII','FontSize',14)
    ylabel('copies ml^{-1}','FontSize',14);set(gca,'FontSize',14)
%     hold on;plot(y,(100/1.5).*y-2200,'r');hold off
    subplot(3,1,3);plot(y,OII(:,t)./(OII(:,t)+OI(:,t)),'k',S(kuro),ostreo(kuro)./100,'ro','LineWidth',2);axis([33.2 34.4 0 1]);title('clade ratio OII/OI','FontSize',14);
    xlabel('salinity','FontSize',14);set(gca,'FontSize',14)
    
    
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
