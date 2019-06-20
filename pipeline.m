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
    [EEGSets, FileInfoSets] = eeglab_clean(fp, files{iFile}, 0, 1, 0);
    
    
    
    
    
end
