function [output, y] = ostreo_model;

% parameters for all experiments
dy = 0.05;
y = 33.2:dy:34.3; % m
k = 5e-08; % diffusity coeff in salinity space, S2/s

% time stepping
dt = 30; % time step in seconds
tend = 10;
tstep = 1:dt:60*60*24*tend;

OI_null=zeros(length(y),length(tstep));
OII_null=zeros(length(y),length(tstep));

OI_pos=zeros(length(y),length(tstep));
OII_pos=zeros(length(y),length(tstep));

OI_neg=zeros(length(y),length(tstep));
OII_neg=zeros(length(y),length(tstep));

% initial conditions and boundary conditions
OI_null(1,:)=1500;
OII_null(end,:)=1500;

OI_pos(1,:)=1500;
OII_pos(end,:)=1500;

OI_neg(1,:)=1500;
OII_neg(end,:)=1500;

% EXPERIMENT 1, unet = 0
u1 = zeros(length(y));
u2 = u1;

m1 = 0;
m2 = m1;

for t = 2:length(tstep);
    for ystep = 2:length(y)-1;
        OI_null(ystep,t) = OI_null(ystep,t-1) + (k/(dy^2)*dt*(OI_null(ystep+1,t-1)-2*OI_null(ystep,t-1)+OI_null(ystep-1,t-1))+u1(ystep)*dt*OI_null(ystep,t-1))-m1*dt*OI_null(ystep,t-1);
        if OI_null(ystep,t)<0, OI_null(ystep,t)=0; end
        OII_null(ystep,t) = OII_null(ystep,t-1) + (k/(dy^2)*dt*(OII_null(ystep+1,t-1)-2*OII_null(ystep,t-1)+OII_null(ystep-1,t-1))+u2(ystep)*dt*OII_null(ystep,t-1))-m2*dt*OII_null(ystep,t-1);
        if OII_null(ystep,t)<0, OII_null(ystep,t)=0; end
    end
end
clear u1 u2 m1 m2

% EXPERIMENT 2, unet > 0
u1 = 0.7/(60*60*24).*ones(length(y));
u2 = u1;

m1 =0.1/(60*60*24);
m2 = m1;

for t = 2:length(tstep);
    for ystep = 2:length(y)-1;
        OI_pos(ystep,t) = OI_pos(ystep,t-1) + (k/(dy^2)*dt*(OI_pos(ystep+1,t-1)-2*OI_pos(ystep,t-1)+OI_pos(ystep-1,t-1))+u1(ystep)*dt*OI_pos(ystep,t-1))-m1*dt*OI_pos(ystep,t-1);
        if OI_pos(ystep,t)<0, OI_pos(ystep,t)=0; end
        OII_pos(ystep,t) = OII_pos(ystep,t-1) + (k/(dy^2)*dt*(OII_pos(ystep+1,t-1)-2*OII_pos(ystep,t-1)+OII_pos(ystep-1,t-1))+u2(ystep)*dt*OII_pos(ystep,t-1))-m2*dt*OII_pos(ystep,t-1);
        if OII_pos(ystep,t)<0, OII_pos(ystep,t)=0; end
    end
end
clear u1 u2 m1 m2

% EXPERIMENT 3, unet < 0
u1 = zeros(length(y));
u2 = u1;

m1 = 0.1/(60*60*24);
m2 = m1;

for t = 2:length(tstep);
    for ystep = 2:length(y)-1;
        OI_neg(ystep,t) = OI_neg(ystep,t-1) + (k/(dy^2)*dt*(OI_neg(ystep+1,t-1)-2*OI_neg(ystep,t-1)+OI_neg(ystep-1,t-1))+u1(ystep)*dt*OI_neg(ystep,t-1))-m1*dt*OI_neg(ystep,t-1);
        if OI_neg(ystep,t)<0, OI_neg(ystep,t)=0; end
        OII_neg(ystep,t) = OII_neg(ystep,t-1) + (k/(dy^2)*dt*(OII_neg(ystep+1,t-1)-2*OII_neg(ystep,t-1)+OII_neg(ystep-1,t-1))+u2(ystep)*dt*OII_neg(ystep,t-1))-m2*dt*OII_neg(ystep,t-1);
        if OII_neg(ystep,t)<0, OII_neg(ystep,t)=0; end
    end
end

output.null = [OI_null(:,t) OII_null(:,t)];
output.pos = [OI_pos(:,t) OII_pos(:,t)];
output.neg = [OI_neg(:,t) OII_neg(:,t)];
