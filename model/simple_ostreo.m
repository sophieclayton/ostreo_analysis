% simple script to calculate concentration of ostreo across front
% Sophie Clayton, October 2011
% sclayton@mit.edu

clear all

% set up domain and set parameters
dy = 1000;
y = 0:dy:30000; % m
k = 0.25*10^4; % m2/s
u1 =(10^-5)-(10^-6); % s-1
u2 =u1;
k1=0.05;
m1=1/(60*60*24*10);
m2=m1;


% time stepping
dt = 30; % quarter day in s
tend = 100;
tstep = 1:dt:60*60*24*tend;

OI=zeros(length(y),length(tstep));
OII=zeros(length(y),length(tstep));

% initial conditions and boundary conditions
OI(1,:)=1000;
OII(end,:)=1000;

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
        OI(ystep,t) = OI(ystep,t-1) + (k/(dy^2)*dt*(OI(ystep+1,t-1)-2*OI(ystep,t-1)+OI(ystep-1,t-1))+u1*dt*OI(ystep,t-1));%-m1*dt*OI(ystep,t-1);
        if OI(ystep,t)<0, OI(ystep,t)=0; end
        OII(ystep,t) = OII(ystep,t-1) + (k/(dy^2)*dt*(OII(ystep+1,t-1)-2*OII(ystep,t-1)+OII(ystep-1,t-1))+u2*dt*OII(ystep,t-1));%-m2*dt*OII(ystep,t-1);
        if OII(ystep,t)<0, OII(ystep,t)=0; end
    end
end

figure;
for t=length(tstep);
    subplot(3,1,1);plot(y,OI(:,t));axis([0 y(end) 0 150]);title(['day =',num2str(tstep(t)/(60*60*24))])
    hold on;plot(y,-100*(1/30000)*y+100,'r');hold off
    subplot(3,1,2);plot(y,OII(:,t));axis([0 y(end) 0 150])
    hold on;plot(y,100*(1/30000)*y,'r');hold off
    subplot(3,1,3);plot(y,OII(:,t)./(OII(:,t)+OI(:,t)).*100);axis([0 y(end) 0 100])
    pause(0.5)
end
