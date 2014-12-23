
function negLL = LLfun(param)
%   this function is called numerous times by fmincon to find the minimum negLL
%   and the associated param vector ([p Ro r g sigma])
%
%   prepared for 7.430 Spring 2008
%   Heidi Sosik, WHOI, Mar. 2008

%%%%%%
%   declare some global variables that will be the same 
%   as in the main script (Riley_optimLL_main)
global I k z1 T Z OneMinusV OneMinusN forcing_time Obs_time Obs_P pc gc

%%%%%%%
%   put your implementation of Riley's model here...
%   this should depend on the vector param ([p Ro r g sigma]) which is passed
%   in to the function each time it is called

time=(0:0.5:375);
irr=interp1(forcing_time,I,time);
oneN=interp1(forcing_time,OneMinusN,time);
oneV=interp1(forcing_time,OneMinusV,time);
temp=interp1(forcing_time,T,time);
zoo=interp1(forcing_time,Z,time);
ksat=interp1(forcing_time,k,time);
zoo1=interp1(forcing_time,z1,time);
time=(0:0.5:375);
p=param(1);
Ro=param(2);
r=param(3);
g=param(4);
sigma=param(5);

dtot=751;
dt=0.5;

phot=NaN(dtot,1);
resp=NaN(dtot,1);
graz=NaN(dtot,1);
dpdt=NaN(dtot,1);
dlnp=NaN(dtot,1);
lnp=NaN(dtot,1);
phy=NaN(dtot,1);
day=NaN(dtot,1);

% initial condition
phy(1)=3.4;
day(1)=0;

% start time stepping
for i=1:dtot;

    phot(i)=((p*irr(i))/(ksat(i)*zoo1(i)))*(1-exp(-ksat(i)*zoo1(i)))*oneN(i)*oneV(i);
    resp(i)=Ro*exp(r*temp(i));
    graz(i)=g*zoo(i);
  
end

for i=2:dtot
    
    dpdt(i)=phy(i-1)*(phot(i-1)-resp(i-1)-graz(i-1));
%    dlnp(i)=0.5*dpdt(i)/phy(i-1);
%    lnp(i)=log(phy(i-1))+dlnp(i);
%    phy(i)=exp(lnp(i));
    phy(i)=phy(i-1)+dpdt(i)*dt;
     
    day(i)=day(i-1)+0.5;
end


%%%%%%
%%%%%%
%   end the function by calculating the negative log-likelihood value (negLL)
%   which depends on both Obs_P and the P estimates from the model run above
%   This is performing a minimization, so make sure you calculate 
%    NEGATIVE log-likeilihood
phy_fit=interp1(sal,phy,Obs_time);
n=7;
sumsq=sum((log(Obs_P)-log(phy_fit)).^2);
negLL=(n/2)*log(2*pi)+(n/2)*log(sigma.^2)+((1/2/sigma.^2)*sumsq);     %<--- your function here
%keyboard
%%%%%
%   use this step to force the optimization to halt if the model blows up
%   thus avoiding an endless search, you'll have to start again with a new
%   initial guess for the parameter vecto
if isnan(negLL)
    negLL = 1e6;
end;
plot(day, phy, Obs_time, Obs_P)