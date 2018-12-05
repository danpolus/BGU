%
%Step 4
%
% https://www.mathworks.com/help/stats/hierarchical-clustering.html
% https://www.mathworks.com/help/stats/dendrogram.html
%optimize by contrast

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select avalanche similarity files','MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

for iFile = 1:length(files)
  AvalancheSimilarityData = load([fp files{iFile}]);

  
  %deal with negative values in similarity matrix
end

save([fp files{iFile}(1:end-4) '_clusters.mat'],'clustering_data');

