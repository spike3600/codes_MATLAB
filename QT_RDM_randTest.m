% the demo script for the function called conditionLabelRandomizationTests4RDMs.m
clear;clc;close all
%% control params
nCond = 12;
nVox = 100;
nSubjects = 15;
options.nRandomisations = 1000;
options.RDMcorrelationType = 'Spearman';%'Kendall_taua';
%% simulate RDMs
refRDM = squareform(pdist(randn(nCond,nVox)));
for subI = 1:nSubjects
    candRDMs(:,:,subI) = squareform(pdist(randn(nCond,nVox)));
end
%% test!
stats_p_r = conditionLabelRandomizationTests4RDMs(refRDM,candRDMs,options)