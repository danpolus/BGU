clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(fullfile(fp, 'runlog.txt'), 'w');
if fid == -1
  error('Cannot open log file.');
end

[files, fp] = uigetfile([fp '*.mat'], 'Select data files', 'MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

fprintf(fid, '%s    pipeline start\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'));
% for iFile = 1:length(files)
parfor iFile = 1:length(files)
    
    fprintf(fid, '%s    1_eeglab_clean: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
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

    fprintf(fid, '%s    2_extract_avalanches: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    AvalancheFileDataSets = extract_avalanches(EEGSets, 0);
    fprintf(fid, '%s    3_get_avalanche_vectors: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    [MultiFileAchVecs, usedTauInfo] = get_avalanche_vectors(AvalancheFileDataSets, 1);
    
%     %load MultiFileAchVecs and usedTauInfo from mat file   
%     [fn, fp] = uigetfile([fp '*.mat'], 'Select avalanche vectors file');
%     load([fp fn],'MultiFileAchVecs','usedTauInfo');

    fprintf(fid, '%s    4_compare_avalanches: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    [SimilarityMat, TestingSet] = compare_avalanches(MultiFileAchVecs, usedTauInfo, 1, 0);
    
%     %load MultiFileAchVecs, SimilarityMat, TestingSet from mat file   
%     [fn, fp] = uigetfile([fp '*.mat'], 'Select similarity matrix file');
%     load([fp fn],'MultiFileAchVecs','SimilarityMat','TestingSet');

    fprintf(fid, '%s    5_cluster_avalanches: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    ClusteringData = cluster_avalanches(MultiFileAchVecs, SimilarityMat, TestingSet, 1, 0);
    
%     %load ClusteringData, MultiFileAchVecs, SimilarityMat, TestingSet from mat file 
%     [fn, fp] = uigetfile([fp '*.mat'], 'Select clustering results file');
%     load([fp fn],'ClusteringData','MultiFileAchVecs','TestingSet');

    fprintf(fid, '%s    6_get_testing_clusters: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    TestingSetClusters = get_testing_clusters(ClusteringData, MultiFileAchVecs, TestingSet, 1);

%     %load TestingSetClusters, ClusteringData from mat file 
%     [fn, fp] = uigetfile([fp '*.mat'], 'Select testing clusters file');
%     load([fp fn],'TestingSetClusters','ClusteringData');

    fprintf(fid, '%s    7_predict_conditions: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    PredictionResults = predict_conditions(TestingSetClusters, ClusteringData, 1, 0);
    
end

fprintf(fid, '%s    pipeline end\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'));
fclose(fid);
