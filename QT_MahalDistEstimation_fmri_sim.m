% QT_MahalDistEstimation_fmri_sim.m
clear;clc;close all
%% settings
nCond = 20;
nVox = 100;
nIter = 1000;
%% compute both estimations for a number of iterations
% Generate a simulationOptions structure.
simulationOptions = simulationOptions_demo_LDt();
[B_true,Y_true,fMRI_a,fMRI_b] = simulateClusteredfMRIData(simulationOptions);
% Let's compute the Mahalanobis distance for one fold, say fold "a" only:
y = fMRI_a.Y;
x = fMRI_a.X;
betas = inv(x'*x)*(x'*y);
res = y - x*betas;
sigma = covdiag(res);
% take two patterns as examples 
pat1 = betas(1,:);
pat2 = betas(2,:);
mahalDist1 = sqrt((pat1-pat2)*inv(sigma)*(pat1-pat2)');
mahalDist2 = pdist([pat1;pat2],'mahalanobis',sigma);
[mahalDist1;mahalDist2]