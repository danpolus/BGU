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
    
%     %load EEGSets from set files
%     [eegSetFiles, fp] = uigetfile([fp '*.set'], 'Select all subject scenario files','MultiSelect','on');
%     if ~iscell(eegSetFiles) %in case only 1 file selected
%         eegSetFiles = {eegSetFiles};
%     end
%     EEGSets = [];
%     for iEegSets = 1:length(eegSetFiles)
%         EEGSets = [EEGSets pop_loadset([fp eegSetFiles{iEegSets}])];
%     end
    
    AvalancheFileDataSets = extract_avalanches(EEGSets, 0);
    [MultiFileAchVecs, usedTauInfo] = get_avalanche_vectors(AvalancheFileDataSets, 1);
    
%     %load MultiFileAchVecs and usedTauInfo from mat file   
%     [fn, fp] = uigetfile([fp '*.mat'], 'Select avalanche vectors file');
%     load([fp fn],'MultiFileAchVecs','usedTauInfo');
    
    [SimilarityMat, TestingSet] = compare_avalanches(MultiFileAchVecs, usedTauInfo, 1, 1);
    
%     %load MultiFileAchVecs, SimilarityMat, TestingSet from mat file   
%     [fn, fp] = uigetfile([fp '*.mat'], 'Select similarity matrix file');
%     load([fp fn],'MultiFileAchVecs','SimilarityMat','TestingSet');

    ClusteringData = cluster_avalanches(MultiFileAchVecs, SimilarityMat, TestingSet, 1, 1);
    
%     %load ClusteringData, MultiFileAchVecs, SimilarityMat, TestingSet from mat file 
%     [fn, fp] = uigetfile([fp '*.mat'], 'Select clustering results file');
%     load([fp fn],'ClusteringData','MultiFileAchVecs','TestingSet');

    TestingSetClusters = get_testing_clusters(ClusteringData, MultiFileAchVecs, TestingSet, 1, 1);

%     %load TestingSetClusters, ClusteringData from mat file 
%     [fn, fp] = uigetfile([fp '*.mat'], 'Select testing clusters file');
%     load([fp fn],'TestingSetClusters','ClusteringData');

    PredictionResults = predict_conditions(TestingSetClusters, ClusteringData, 1, 1);
    
end
