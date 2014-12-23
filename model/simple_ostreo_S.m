% simple script to calculate concentration of ostreo across front using a
% simple reaction-diffusion model
% Sophie Clayton, October 2011, updated December 2014
% sclayton@mit.edu, sclayton@uw.edu

%clear all

n=24;

load ../data/ostreo % load the ostreo clade abundance data

kuro=find(lon>140 & z==0); % use only the data from the Kuroshio

[out, y] = ostreo_model;

run('../figures/plot_ostreo_model.m')