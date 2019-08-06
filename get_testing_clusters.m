%
%Step 7 : get_testing_clusters - find clusters for training and testing set vectors
%
%  inputs:
% ClusteringDataSets - clusters and statistics
% MultiFileAchVecs - avalanche vectors from all sets
% SimilarityMat - similarity between avalanches distances matrix
% TrainingSets - avalanches (epochs) for training
% TestingSets - avalanches (epochs) for testing
% save_str - add this sting to filename
% saveFlg
%
%  outputs:
% TrainingClusterSets - training set clusters
% TestingClusterSets- testing set clusters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TrainingClusterSets, TestingClusterSets] = get_testing_clusters(ClusteringDataSets, MultiFileAchVecs, SimilarityMat, TrainingSets, TestingSets, save_str, saveFlg)

optimal_tau_t = MultiFileAchVecs{1}(~cellfun(@isempty,{MultiFileAchVecs{1}.is_optimal_tau}));
fileInfo = optimal_tau_t(([MultiFileAchVecs{1}.is_optimal_tau] == 1)).dataInfo.FileInfo;

if saveFlg
    output_fp = [fileInfo.base_fp '5 testing\'];
    mkdir(output_fp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TrainingClusterSets = [];
TestingClusterSets = [];

percent_waitbar = 0;
f_waitbar = waitbar(percent_waitbar, ['get testing set clusters ' num2str(100*percent_waitbar) '%'], 'Name',fileInfo.orig_fn );
                
for iCrossValid = 1:length(TrainingSets)
    for iTau = 1:length(ClusteringDataSets{iCrossValid})
        if isempty(ClusteringDataSets{iCrossValid}(iTau).tau)
            continue;
        end
        
        percent_waitbar = iCrossValid/length(TrainingSets);
        waitbar(percent_waitbar,f_waitbar,['cross-valid ' num2str(iCrossValid) '  get testing set clusters  ' num2str(100*percent_waitbar) '%']);
            
        TrainingClusterSets{iCrossValid}(iTau) = prepare_clusters(MultiFileAchVecs, SimilarityMat, ClusteringDataSets{iCrossValid}, TrainingSets{iCrossValid}, fileInfo, iTau, 1);
        TestingClusterSets{iCrossValid}(iTau) = prepare_clusters(MultiFileAchVecs, SimilarityMat, ClusteringDataSets{iCrossValid}, TestingSets{iCrossValid}, fileInfo, iTau, 0);       
    end
end

close(f_waitbar);

if saveFlg
    save([output_fp fileInfo.orig_fn '_' save_str 'Testclust.mat'],'TestingClusterSets','TrainingClusterSets','ClusteringDataSets');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetClusters = prepare_clusters(MultiFileAchVecs, SimilarityMat, ClusteringData, TestingSet, fileInfo, iTau, isTrain)

%testing set
SetClusters.CondClst = [];
SetClusters.CondIds = TestingSet(iTau).CondIds;
SetClusters.tau = TestingSet(iTau).tau;
SetClusters.is_optimal_tau = TestingSet(iTau).is_optimal_tau;
SetClusters.fileInfo.base_fp = fileInfo.base_fp;
SetClusters.fileInfo.orig_fn = fileInfo.orig_fn;

for iCond = 1:length(TestingSet(iTau).CondIds)
    for iEpoch = 1:length(TestingSet(iTau).EpochIds{iCond})
        fln = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(4:6));
        epc = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(10:12));
        ach_vectors_t = MultiFileAchVecs{fln}(iTau).epochs_vecs{epc};
        for iVec = 1:length(ach_vectors_t)
            if isTrain
                [SetClusters.CondClst(iCond).EpochClst(iEpoch).cluster_num(iVec,:), ...
                    SetClusters.CondClst(iCond).EpochClst(iEpoch).cluster_sim(iVec,:)] = ...
                    find_training_avalanche_cluster(ach_vectors_t(iVec), ClusteringData(iTau));
            else
                [SetClusters.CondClst(iCond).EpochClst(iEpoch).cluster_num(iVec,:), ...
                    SetClusters.CondClst(iCond).EpochClst(iEpoch).cluster_sim(iVec,:)] = ...
                    find_testing_avalanche_cluster(ach_vectors_t(iVec), ClusteringData(iTau), SimilarityMat(iTau));
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cluster_num, cluster_sim] = find_training_avalanche_cluster(avalanche_vec, ClusteringData)

nofLen = length(ClusteringData.Clusters);
cluster_num = zeros(1,nofLen);
cluster_sim = zeros(1,nofLen);
ach_idx = contains(ClusteringData.Id{avalanche_vec.length_bins}, avalanche_vec.id);
cluster_num(avalanche_vec.length_bins) = ClusteringData.Clusters{avalanche_vec.length_bins}(ach_idx);
ach_idx = contains(ClusteringData.Id{nofLen-1}, avalanche_vec.id);
cluster_num(nofLen-1) = ClusteringData.Clusters{nofLen-1}(ach_idx);
ach_idx = contains(ClusteringData.Id{nofLen}, avalanche_vec.id);
cluster_num(nofLen) = ClusteringData.Clusters{nofLen}(ach_idx);
cluster_sim(cluster_num > 0) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cluster_num, cluster_sim] = find_testing_avalanche_cluster(avalanche_vec, ClusteringData, SimilarityMat)

nofLen = length(ClusteringData.Clusters);
cluster_num = zeros(1,nofLen);
cluster_sim = zeros(1,nofLen);

lenInx = [avalanche_vec.length_bins  nofLen-1  nofLen];
for iLen = lenInx
    
    clusters = unique(ClusteringData.Clusters{iLen});
    if isempty(clusters)
        continue;
    end
    
    max_sim = zeros(1,length(clusters));
    med_sim = zeros(1,length(clusters));
    mahalanobis_med_sim = zeros(1,length(clusters));
    for iClust = 1:length(clusters)
        ClusterAchIds = ClusteringData.Id{iLen}(clusters(iClust) == ClusteringData.Clusters{iLen});
        if iLen == nofLen %concat
            iLenMat = avalanche_vec.length_bins;
            if isempty(ClusteringData.Id{iLenMat})
                continue;
            end               
            ClusterAchIds = ClusterAchIds(contains(ClusterAchIds, ClusteringData.Id{iLenMat}));
            if isempty(ClusterAchIds)
                continue;
            end            
        else
            iLenMat = iLen;
        end
        ClusterIds_idx = contains(SimilarityMat.Id{iLenMat}, ClusterAchIds);
        ach_idx = contains(SimilarityMat.Id{iLenMat}, avalanche_vec.id);
        similarities = SimilarityMat.Mat{iLenMat}(ClusterIds_idx,ach_idx);
        
        max_sim(iClust) = max(similarities);
        if max_sim(iClust) == 1 %exact match -> cluster found
            med_sim(iClust) = 1;
            mahalanobis_med_sim(iClust) = Inf;
            break;
        else
            med_sim(iClust) = max(median(similarities),mean(similarities));
            std_sim = std(similarities);
            if std_sim > 0
                mahalanobis_med_sim(iClust) = med_sim(iClust)/std_sim;
            end
        end
    end
    [~,max_inx] = max(mahalanobis_med_sim);%max(med_sim);max(max_sim);
    cluster_num(iLen) = clusters(max_inx);
    cluster_sim(iLen) = med_sim(max_inx); %max_sim(max_inx);

end
