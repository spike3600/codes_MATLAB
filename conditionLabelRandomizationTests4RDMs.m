function stats_p_r = conditionLabelRandomizationTests4RDMs(refRDM,candRDMs,options)
% compute the correlation of the refRDM to candRDMs and test the
% significance of the correlations via condition-label randomization test.
% Output is a structure containing the correlations and the p-values 
% (uncorrected and also FWE corrected).
% the function performs a condition-label randomisation test on the
% relatedness of the refRDM to all of the RDMs in the candRDMs.
% INPUTs:
% refRDM: an [n n] RDM that will be tested against all of the candidates
% candRDMs: [n n nCandidates] matrix
% options is a strcture containing two fields:
%           options.nRandomisations: number of permutations applied
%           options.RDMcorrelationType: the measure of the RDM similarity.
%           The code also handles Kendall_taua which is the prefered
%           measure when there are tied ranks in either of the RDMs.

% Hamed Nili
% note that due to code dependencies, you'll need to add the path to the 
% rsatoolbox before running the code. 

nCandRDMs = size(candRDMs,3);
nRandomisations = options.nRandomisations;

for candI = 1:nCandRDMs
        if isequal(options.RDMcorrelationType,'Kendall_taua')
            cand2refSims(candI)=rankCorr_Kendall_taua(vectorizeRDMs(candRDMs(:,:,candI))',vectorizeRDMs(refRDM)');
        else
            cand2refSims(candI)=corr(vectorizeRDMs(candRDMs(:,:,candI))',vectorizeRDMs(refRDM)','type',options.RDMcorrelationType,'rows','pairwise');
        end
end

% vectorize all the candidate RDMs
for rdmI = 1:nCandRDMs
    rdms(rdmI,:) = vectorizeRDM(candRDMs(:,:,rdmI));
end
[n,n]=size(refRDM);
exhaustPermutations = false;
if n < 8
    allPermutations = exhaustivePermutations(n);
    nRandomisations = size(allPermutations, 1);
    exhaustPermutations = true;
    warning('(!) Comparing RDMs with fewer than 8 conditions (per conditions set) will produce unrealiable results!\n  + I''ll partially compensate by using exhaustive instead of random permutations...');
end%if n < 8
% make space for null-distribution of correlations
rs_null=nan(options.nRandomisations,nCandRDMs);

%tic
if isequal(options.RDMcorrelationType,'Kendall_taua')
    for randomisationI=1:options.nRandomisations
        if exhaustPermutations
            randomIndexSeq = allPermutations(randomisationI, :);
        else
            randomIndexSeq = randomPermutation(n);
        end%if
        
        rdmA_rand_vec=vectorizeRDM(refRDM(randomIndexSeq,randomIndexSeq));
        for candI = 1:nCandRDMs
            rs_null(randomisationI,candI)=rankCorr_Kendall_taua(rdmA_rand_vec',rdms(candI,:)');
        end
        if mod(randomisationI,floor(options.nRandomisations/100))==0
            fprintf('%d%% ',floor(100*randomisationI/options.nRandomisations))
            if mod(randomisationI,floor(options.nRandomisations/10))==0, fprintf('\n'); end;
        end
    end % randomisationI
    fprintf('\n');
else
    for randomisationI=1:options.nRandomisations
        if exhaustPermutations
            randomIndexSeq = allPermutations(randomisationI, :);
        else
            randomIndexSeq = randomPermutation(n);
        end%if
        
        rdmA_rand_vec=vectorizeRDM(refRDM(randomIndexSeq,randomIndexSeq));
        rs_null(randomisationI,:)=corr(rdmA_rand_vec',rdms','type',options.RDMcorrelationType,'rows','pairwise');
        if mod(randomisationI,floor(options.nRandomisations/100))==0
            fprintf('%d%% ',floor(100*randomisationI/options.nRandomisations))
            if mod(randomisationI,floor(options.nRandomisations/10))==0, fprintf('\n'); end;
        end
    end % randomisationI
    fprintf('\n');
end

% p-values from the randomisation test
for candI = 1:nCandRDMs
    p_randCondLabels(candI) = 1 - relRankIn_includeValue_lowerBound(rs_null(:,candI),cand2refSims(candI)); % conservative
end
% p-values corrected for MC by controlling FWE rate
for candI = 1:nCandRDMs
    p_randCondLabels_fwe(candI) = 1 - relRankIn_includeValue_lowerBound(max(rs_null'),cand2refSims(candI)); % conservative
end
stats_p_r.candRelatedness_p_uncorr = p_randCondLabels;
stats_p_r.candRelatedness_p_fwe = p_randCondLabels_fwe;
stats_p_r.candRelatedness_r = cand2refSims;
stats_p_r.nullRs_r = rs_null;