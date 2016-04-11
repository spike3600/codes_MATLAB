function [] = writeMaps(map,binaryMask,structPath,resultsPath,readFile,options)
% this is a stand-alone function that saves a map (e.g. statistical map) in
% various SPM-compatible formats. INPUTS: map: 3D map (e.g. t-statistic of
% a contrast or the correlation ... map from SL analysis). binaryMask: a
% binary mask that would be applied to the input map. structPath: path to
% the folder containing structural data of the subject for which the map is
% computed. resultsPath: path to where the resulting maps and masks would
% be saved options: a structure specifying the operation(s) done by the
% function
%   options.voxelSize: a triple indicating the size of the voxels in mm
%   options.analysisName: a string containing the name of the analysis
%   performed. Subject name, model name, etc. will be specified here.
%   options.writeNative: a bolean variable. If set to 1, the native space
%   maps and mask would be saved as image files. options.writeNormalised: a
%   bolean variable. If set to 1, the native space maps and mask would be
%   warped and saved as image files. options.writeSmoothedNormalised: a
%   bolean variable. If set to 1, the native space maps and mask would be
%   warped,smoothed, masked and then saved as image files.

% Hamed Nili, April 2016 (inspired by the code written by Cai Wingfield)
%% set options
options = setIfUnset(options,'voxelSize',[2 2 2]);
options = setIfUnset(options,'analysisName','unnamed');
options = setIfUnset(options,'writeNative',1);
options = setIfUnset(options,'writeNormalised',1);
options = setIfUnset(options,'writeSmoothedNormalised',1);
options = setIfUnset(options,'fwhm',[5 5 5]);
if options.writeSmoothedNormalised, options.writeNormalised = 1;end
if nargin == 2,options = struct();end % if: nargin
if ~exist('map','var')||numel(size(map)) ~= 3,error('you need to pass a 3D volume as the input map');end
if ~exist('structPath','var')&&(options.writeNormalised == 1 || options.writeSmoothedNormalised == 1),
    error('you need to provide the path to structural images');end
if ~exist('readFile','var')&&(options.writeNormalised == 1 || options.writeSmoothedNormalised == 1),
    error('you need to provide the path to a pre-processed SPM file of the same subject');end
if ~exist('resultsPath','var'),error('you need to provide the results path');end
if ~isequal(size(map),size(binaryMask)),
    error('the first two inputs need to have the same size');
end
warpFlags.interp = 1;
warpFlags.wrap = [0 0 0];
warpFlags.vox = options.voxelSize; 
warpFlags.bb = [-78 -112 -50; 78 76 85];
warpFlags.preserve = 0;

%% write the native-space map and binary mask
cd(structPath);
% load the header information
V = spm_vol(readFile);

% write the native-space masked map to disk
mapMetadataStruct_nS = V;
mapMetadataStruct_nS.fname = fullfile(resultsPath,[options.analysisName,'_nativeSpaceMap.img']);
mapMetadataStruct_nS.descrip =  'map';
mapMetadataStruct_nS.dim = size(map);
maskMetadataStruct_nS = V;

if options.writeNative
    
    fprintf('writing the native space map...\n')
    spm_write_vol(mapMetadataStruct_nS, map);
    fprintf('... done\n')
    
    % write the native space binary mask to disk
    maskMetadataStruct_nS.fname = fullfile(resultsPath,[options.analysisName,'_nativeSpaceBinaryMask.img']);
    maskMetadataStruct_nS.descrip =  'Native space mask';
    maskMetadataStruct_nS.dim = size(binaryMask);
    
    fprintf('writing the native space binary mask...\n')
    spm_write_vol(maskMetadataStruct_nS, binaryMask);
    fprintf('... done\n')
end
if options.writeNormalised
    % Load in common space warp definition
    warpDefFilenames = dir(fullfile(structPath,'*_seg_sn.mat'));
    warpDefFilename = warpDefFilenames(end).name;
    % Warp and write common space maps to disk
    fprintf('writing the normalised map...\n')
    spm_write_sn(mapMetadataStruct_nS,warpDefFilename,warpFlags);% write the normalised map to disk
    fprintf('... done\n')
    
    % Warp and write common space masks to disk
    fprintf('writing the normalised binary mask...\n')
    spm_write_sn(maskMetadataStruct_nS,warpDefFilename,warpFlags);% write the normalised mask to disk
    fprintf('... done\n')
    
    % Now read them back in
    
    % Where are they?
    [warpedPath_map, warpedFile_map, warpedExt_map] = fileparts(mapMetadataStruct_nS.fname);
    [warpedPath_mask, warpedFile_mask, warpedExt_mask] = fileparts(maskMetadataStruct_nS.fname);
    
    % Warped versions are prefixed with 'w'
    warpedFile_map = ['w' warpedFile_map];
    warpedFile_mask = ['w' warpedFile_mask];
    
    % Read the normalised mask from the disk
    mask_sS = spm_read_vols(spm_vol(fullfile(warpedPath_mask, [warpedFile_mask warpedExt_mask])));
    
    % Fix the normalisation of the mask
    maskMetadataStruct_sS = spm_vol(fullfile(warpedPath_map, [warpedFile_map warpedExt_map]));
    %     maskMetadataStruct_sS.fname = fullfile(userOptions.rootPath, 'Maps', [userOptions.analysisName '_commonSpaceMask_' maskName '_' modelName '_' subject '.img']);
    maskMetadataStruct_sS.fname = fullfile(resultsPath,[options.analysisName '_commonSpaceMask.img']);
    maskMetadataStruct_sS.descrip =  'Common space mask';
    maskMetadataStruct_sS.dim = size(mask_sS);
    
    maskThreshold = 0.01;
    mask_sS(mask_sS < maskThreshold) = 0;
    mask_sS(isnan(mask_sS)) = 0;% note that this is not *perfect* and
    % generally there might be other ways of tackling the question which
    % obviate the need to normalise binary masks.
    
    maskMetadataStruct_sS.dim = size(mask_sS);
    fprintf('writing the ''fixed'' normalised map, called the common space mask...\n')
    spm_write_vol(maskMetadataStruct_sS, mask_sS);
    fprintf('... done\n')
end

if options.writeSmoothedNormalised
    % Smooth the normalised map
    
    [warpedPath_map, warpedFile_map, warpedExt_map] = fileparts(mapMetadataStruct_nS.fname);
    [warpedPath_mask, warpedFile_mask, warpedExt_mask] = fileparts(maskMetadataStruct_nS.fname);
    
    % Warped versions are prefixed with 'w'
    warpedFile_map = ['w' warpedFile_map];
    warpedFile_mask = ['w' warpedFile_mask];
    
    % Smoothed versions are prefixed with 's'
    smoothedWarpedFile_map = ['s' warpedFile_map];
    
    % Smooth it
    smoothingKernel_fwhm = options.fwhm;
    spm_smooth(fullfile(warpedPath_map, [warpedFile_map warpedExt_map]), fullfile(warpedPath_map, [smoothedWarpedFile_map warpedExt_map]), smoothingKernel_fwhm);
    
    % Read it back in
    smoothedDataMetadataStruct = spm_vol(fullfile(warpedPath_map, [smoothedWarpedFile_map warpedExt_map]));
    smoothedData = spm_read_vols(smoothedDataMetadataStruct);
    
    % Mask the smoothed data by the sS mask
    maskedData = smoothedData;
    maskedData(mask_sS == 0) = NaN;
    %     maskedSmoothedRMaps_sS.(modelName).(subject).(maskName) = maskedData;
    
    % Write it back to disk
    maskedDataMetadataStruct_sS = smoothedDataMetadataStruct;
    maskedDataMetadataStruct_sS.fname = fullfile(resultsPath,['msw' options.analysisName '.img']); % 'msw' for 'masked, smoothed, warped'
    maskedDataMetadataStruct_sS.descrip =  'Masked smoothed normalised data';
    maskedDataMetadataStruct_sS.dim = size(maskedData);
    
    fprintf('writing the masked smoothed normalised map (msw prefix)...\n')
    spm_write_vol(maskedDataMetadataStruct_sS, maskedData);
    fprintf('... done\n')
    
end

cd(resultsPath);










