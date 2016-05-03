% QT_parCorrVsRegression.m
% When comparing a reference RDM to a candidate RDM, we might want to
% partial out the contribution of another RDM on the correlation between
% the candidate and the reference RDMs.
% this could be done via partial correlations or linear regression.
% this scripts demonstrates that the two approaches are identical.
% 
clear;clc;close all
%% control parameters
nCond = 72;
nVox = 100;

%% simulating the RDMs
rdm_ref = pdist(randn(nCond,nVox))';
rdm_cand1 = pdist(randn(nCond,nVox))';
rdm_cand2 = pdist(randn(nCond,nVox))';
%% zscore all the RDMs
rdm_ref = zscore(rdm_ref);
rdm_cand1 = zscore(rdm_cand1);
rdm_cand2 = zscore(rdm_cand2);
%% partial correlation
[r1,p1] = partialcorr(rdm_ref,rdm_cand1,rdm_cand2);
fprintf('partial correlation coefficient = %.3f \n',r1)
%% linear regression
computeRes = @(x,y) y-(x*(x\y));% residuals from GLM
% regress-out rdm2 from the reference rdm
res_ref = computeRes(rdm_cand2,rdm_ref);
% regress-out rdm2 from the candidate rdm
res_cand = computeRes(rdm_cand2,rdm_cand1);

[r2,p2] = corr(res_cand,res_ref);% this is equal to res_ref\res_cand
fprintf('regression-based results = %.3f \n',r2)


