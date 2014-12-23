% Riley_optimLL_main.m runs an optimization to find estimates of model parameters of
% given the set of observations reported in Riley (1946).
% This scripts calls fmincon, which in turns call the user-defined LLfun.m
% LLfun.m must be edited to specify the plankton model and associated
% log-likelihood function.
% prepared for 7.430 Spring 2008
% Heidi Sosik, WHOI, Mar. 2008

% declare some global variables that will also appear inside LLfun
global I k z1 T Z OneMinusV OneMinusN forcing_time Obs_time Obs_P 

warning off

% load the forcing data and the observations
load ostreo


%%%%%%
%   set initial guess for the parameter vector; 
%   Riley's model blows up easily so a limited range on some parameters 
%   can help to get started in a valid solution
% param_init = [3 0.0175  0.069 0.0075 1].*rand(5,1)';  % p Ro r g sigma

param_init = [2.5 .0175 .069 .0075 0.1591]

%%%%%%
% lower constraints on parameters
LB = [0 0 0 0 0]; %p Ro r g sigma
% upper constraints on parameters, 
%   upper limits help to speed optimization and avoid model blow up
UB = [10 5 5 5 5]; % p Ro r g sigma

%%%%%%
% set some options for the optimization
options = optimset('MaxIter', 100000, 'Display', 'iter', 'MaxFunEvals', 100000, 'TolFun', 1e-8, 'TolX', 1e-8, 'TolCon', 1e-12);

%%%%%
%   call the optimization function
%   LLfun is your custom function (separate m-file) that runs the model and
%   calculates negative log-likelihood, given parameter values as input
[param_estimated,negLL_final] = fmincon('LLfun', param_init,[],[],[],[],LB,UB,[],options);
%   param_estimated ([p Ro r g sigma]) are the parameter values that provide
%   the minimum negative log-likelihood
%   IF the optimization terminates in an acceptable manner - Be careful about
%   what you get!

% Use this as a simple-minded way to start again (with another initial guess) if the model blows up
if negLL_final == 1e6,  %presumes case where LLfun resets negLL = NaN (model blow up) to negLL = 1e6
    ostreo_optimLL_main
end;
