clear;clc
%% control params
nCond = 12;
nVox = 100;
nSubjects = 15;
options.nRandomisations = 1000;
options.RDMcorrelationType = 'Kendall_taua';
%% simulate RDMs
refRDM = squareform(pdist(randn(nCond,nVox)));
for subI = 1:nSubjects
    candRDMs(:,:,subI) = squareform(pdist(randn(nCond,nVox)));
end
%% test!
stats_p_r = conditionLabelRandomizationTests4RDMs(refRDM,candRDMs,options)


