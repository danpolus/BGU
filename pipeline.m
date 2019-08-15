clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

saveLog = 1;
profile off

%when need to load 2 files from different folders
middleRun = [];
middleRun.sufix1 = '_avalanches.mat';
middleRun.sufix2 = '_similarity.mat';
middleRun.folder2 = '3 similarity';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveLog
    fid = fopen(fullfile(fp, 'runlog.txt'), 'w');
    if fid == -1
        error('Cannot open log file.');
    end
end

[files, fp] = uigetfile([fp '*.mat'], 'Select data files', 'MultiSelect','on');
if ~iscell(files) %in case only 1 file selected
    files = {files};
end

if ~isempty(middleRun)
    k = strfind(fp,'\');
    middleRun.base_fp = fp(1:k(end-1));
end

if saveLog
    fprintf(fid, '%s    pipeline start\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'));
end
for iFile = 1:length(files)
% parfor iFile = 1:length(files)
    
    if isempty(middleRun) %run full pipeline
        
        if saveLog
            fprintf(fid, '%s    1_eeglab_clean: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
        end
        EEGSets = eeglab_clean(fp, files{iFile}, 0, 1, 0);
        
%         %load EEGSets from set files
%         [eegSetFiles, fp] = uigetfile([fp '*.set'], 'Select all subject scenario files','MultiSelect','on');
%         if ~iscell(eegSetFiles) %in case only 1 file selected
%             eegSetFiles = {eegSetFiles};
%         end
%         EEGSets = [];
%         for iEegSets = 1:length(eegSetFiles)
%             EEGSets = [EEGSets pop_loadset([fp eegSetFiles{iEegSets}])];
%         end
        
        if saveLog
            fprintf(fid, '%s    2_extract_avalanches: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
        end
        AvalancheFileDataSets = extract_avalanches(EEGSets, 0);
        if saveLog
            fprintf(fid, '%s    3_get_avalanche_vectors: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
        end
        [MultiFileAchVecs, usedTauInfo] = get_avalanche_vectors(AvalancheFileDataSets, 1);
%         MultiFileAchVecs = simulate_avalanche_vectors(MultiFileAchVecs, usedTauInfo, 1); %simulate artificial input
        
%         %load MultiFileAchVecs and usedTauInfo from mat file
%         [fn, fp] = uigetfile([fp '*.mat'], 'Select avalanche vectors file');
%         load([fp fn],'MultiFileAchVecs','usedTauInfo');
        
    if saveLog
        fprintf(fid, '%s    4_compare_avalanches: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    end
    SimilarityMat = compare_avalanches(MultiFileAchVecs, usedTauInfo, 1, 0);
    
%     %load MultiFileAchVecs, SimilarityMat from mat file
%     [fn, fp] = uigetfile([fp '*.mat'], 'Select similarity matrix file');
%     load([fp fn],'MultiFileAchVecs','SimilarityMat');


    else %if isempty(middleRun)        
        load([fp files{iFile}],'MultiFileAchVecs','usedTauInfo');
        base_fn = files{iFile}(1:end-length(middleRun.sufix1));
        load([middleRun.base_fp  middleRun.folder2 '\' base_fn  middleRun.sufix2],'MultiFileAchVecs','SimilarityMat');
    end
    
    if saveLog
        fprintf(fid, '%s    5_create_train_test_sets: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    end
    TrainValidTest = create_train_test_sets(MultiFileAchVecs, usedTauInfo);
    if saveLog
        fprintf(fid, '%s    6_cluster_avalanches: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    end
    ValidationClusteringDataSets = cluster_avalanches(MultiFileAchVecs, SimilarityMat, {TrainValidTest.CrossValid.TrainingSet}, TrainValidTest, 'validTrain', 1, 0);
    FinalClusteringDataSet = cluster_avalanches(MultiFileAchVecs, SimilarityMat, {TrainValidTest.TrainingSet}, TrainValidTest, 'finalTrain', 1, 0);

%     %load ClusteringDataSets, MultiFileAchVecs, TrainValidTest from mat file
%     [fn, fp] = uigetfile([fp '*.mat'], 'Select clustering results file');
%     load([fp fn],'ClusteringDataSets','MultiFileAchVecs','TrainValidTest');
%     ValidationClusteringDataSets = ClusteringDataSets;
%     FinalClusteringDataSet = ClusteringDataSets;
    
    if saveLog
        fprintf(fid, '%s    7_get_testing_clusters: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    end
    [ValidationTrainingClusterSets, ValidationTestingClusterSets] = get_testing_clusters(ValidationClusteringDataSets, MultiFileAchVecs, SimilarityMat, {TrainValidTest.CrossValid.TrainingSet}, {TrainValidTest.CrossValid.TestingSet});
    [FinalTrainingClusterSet, FinalTestingClusterSet] = get_testing_clusters(FinalClusteringDataSet, MultiFileAchVecs, SimilarityMat, {TrainValidTest.TrainingSet}, {TrainValidTest.TestingSet});
    if saveLog
        fprintf(fid, '%s    8_predict_conditions: %s\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'), files{iFile});
    end
    ValidTrainPredictionResults = predict_conditions(ValidationTrainingClusterSets, ValidationClusteringDataSets, 'validTrain', 1, 0);
    ValidTestPredictionResults = predict_conditions(ValidationTestingClusterSets, ValidationClusteringDataSets, 'validTest', 1, 0);
    FinalTrainPredictionResults = predict_conditions(FinalTrainingClusterSet, FinalClusteringDataSet, 'finalTrain', 1, 0);
    FinalTestPredictionResults = predict_conditions(FinalTestingClusterSet, FinalClusteringDataSet, 'finalTest', 1, 0);
    
end

if saveLog
    fprintf(fid, '%s    pipeline end\n', datestr(now, 'yyyy/mm/dd HH:MM:SS.FFF'));
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile viewer
profile off
