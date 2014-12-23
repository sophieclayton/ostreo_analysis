function [output] = ostreo_model(par);

% time stepping
dt = 30; % time step in seconds
tend = 10;
tstep = 1:dt:60*60*24*tend;

OI=zeros(length(y),length(tstep));
OII=zeros(length(y),length(tstep));

% initial conditions and boundary conditions
OI(1,:)=5000;
OII(end,:)=5000;

% OI(1,2:end)=0.1;
% OII(end,2:end)=0.1;

% set gaussian nutrient profile
%k=10^4n0=5;
% n=0.5 + n0*exp(-((y-15000).^2)./(2*5000^2));
% umax = 1/(60*60*24);
% u1=(umax.*(n./(n+k1)));
% u2=u1;
%n=ones(length(y),length(tstep));

for t = 2:length(tstep);
%    n0=5*(1+sin((2*pi*1/(60*60*24*3))*tstep(t)));
%    n(:,t)=0.5 + n0*exp(-((y-15000).^2)./(2*2000^2));
%    u1=(umax.*(n(:,t)./(n(:,t)+k1)));
%    u2=u1;
    for ystep = 2:length(y)-1;
        OI(ystep,t) = OI(ystep,t-1) + (k/(dy^2)*dt*(OI(ystep+1,t-1)-2*OI(ystep,t-1)+OI(ystep-1,t-1))+u1(ystep)*dt*OI(ystep,t-1))-m1*dt*OI(ystep,t-1);
        if OI(ystep,t)<0, OI(ystep,t)=0; end
        OII(ystep,t) = OII(ystep,t-1) + (k/(dy^2)*dt*(OII(ystep+1,t-1)-2*OII(ystep,t-1)+OII(ystep-1,t-1))+u2(ystep)*dt*OII(ystep,t-1))-m2*dt*OII(ystep,t-1);
        if OII(ystep,t)<0, OII(ystep,t)=0; end
    end
end

output = [OI(:,tend) OII(:,tend)];
