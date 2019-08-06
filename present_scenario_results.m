clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[files, fp] = uigetfile([fp '*.mat'], 'Select results data files', 'MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

allSubjectsResults = [];
for iFile = 1:length(files)
    
    load([fp files{iFile}],'PredictionResultSets','ClusteringDataSets');
    
    tp_total = zeros(4);
    Kappa = 0;
    for iCrossValid = 1:length(PredictionResultSets)
        Kappa = Kappa + 1 - (1-roc.tp_total)/(1-Stats.P_cond);
        tp_total = tp_total + roc.tp_total;
    end
    allSubjectsResults(iFile).Kappa = Kappa/length(PredictionResultSets);
    allSubjectsResults(iFile).tp_total = tp_total/length(PredictionResultSets); 
    
end

plot([allSubjectsResults.Kappa]);
