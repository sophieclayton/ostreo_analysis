% simple script to calculate concentration of ostreo across front using a
% simple reaction-diffusion model
% Sophie Clayton, October 2011, updated December 2014
% sclayton@mit.edu, sclayton@uw.edu

%clear all

n=24;

load ../data/ostreo % load the ostreo clade abundance data

kuro=find(lon> 140 & z==0); % use only the data from the Kuroshio

[out, y] = ostreo_model;

run('../figures/plot_ostreo_model.m')

% estimate fit to observations
s_range=find(S>=33.2 & S<=34.3);
S_int = sort(S(s_range));
obs = O(s_range,:);
obs = [obs(:,1); obs(:,2)];

mod_pos = [interp1(y,out.pos(:,1),S_int); interp1(y,out.pos(:,2),S_int)];
mod_neg = [interp1(y,out.neg(:,1),S_int); interp1(y,out.neg(:,2),S_int)];
mod_null  = [interp1(y,out.null(:,1),S_int); interp1(y,out.null(:,2),S_int)];

error = zeros(3,1);
error(1) = sqrt(mean((obs-mod_pos).^2));
error(2) = sqrt(mean((obs-mod_neg).^2));
error(3) = sqrt(mean((obs-mod_null).^2));

error

