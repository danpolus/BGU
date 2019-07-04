clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select data files', 'MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

% for iFile = 1:length(files)
parfor iFile = 1:length(files)
    EEGSets = eeglab_clean(fp, files{iFile}, 0, 1, 0);
    
%     [eegSetFiles, fp] = uigetfile([fp '*.set'], 'Select all subject scenario files','MultiSelect','on');
%     if ~iscell(eegSetFiles) %in case only 1 file selected
%         eegSetFiles = {eegSetFiles};
%     end
%     EEGSets = [];
%     for iEegSets = 1:length(eegSetFiles)
%         EEGSets(iEegSets) = pop_loadset([fp eegSetFiles{iEegSets}]);
%     end
    
    AvalancheFileDataSets = extract_avalanches(EEGSets, 0);
    
    
    
    
end
