% QT_writeMaps.m
clear;clc;close all hidden
%% load a sample t-map
exampleDir = '/Volumes/SPIKE/Oxford/fMRI/blah/structurals/del';
exampleFile = fullfile(exampleDir,'spmT_0013.nii');
readFile = exampleFile;
map = spm_read_vols(spm_vol(exampleFile));
maskFile = fullfile(exampleDir,'mask.nii');
binaryMask = spm_read_vols(spm_vol(maskFile));
structPath = '/Volumes/SPIKE/Oxford/fMRI/blah/structurals';
options.analysisName = 'BLAH';
options.writeNative = 1;
resultsPath = '/Volumes/SPIKE/Oxford/fMRI/blah/structurals/del';

% create a template header info
% V.fname = '';
% V.mat = eye(4);
% V.dim = [0 0 0 0];
% V.dt = [16 0];
% V.pinfo = [1;0;0];% this is a spawn SPM structure
% V = spm_vol(exampleFile)
% V.pinfo = [1;0;0];
% V.mat = eye(4);
% write the native-space masked map to disk
% mapMetadataStruct_nS = V;
% mapMetadataStruct_nS.fname = fullfile(resultsPath,[options.analysisName,'_nativeSpaceMap.img']);
% mapMetadataStruct_nS.descrip =  'map';
% mapMetadataStruct_nS.dim = size(map);
% maskMetadataStruct_nS = V;
% 
% if options.writeNative
%     
%     fprintf('writing the native space map...\n')
%     spm_write_vol(mapMetadataStruct_nS, map);
%     fprintf('... done\n')
% end

writeMaps(map,binaryMask,structPath,resultsPath,readFile,options)