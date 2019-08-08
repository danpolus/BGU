%
%Step 7 : get_testing_clusters - find clusters for training and testing set vectors
%
%  inputs:
% ClusteringDataSets - clusters and statistics
% MultiFileAchVecs - avalanche vectors from all sets
% SimilarityMat - similarity between avalanches distances matrix
% TrainingSets - avalanches (epochs) for training
% TestingSets - avalanches (epochs) for testing
%
%  outputs:
% TrainingClusterSets - training set clusters
% TestingClusterSets- testing set clusters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TrainingClusterSets, TestingClusterSets] = get_testing_clusters(ClusteringDataSets, MultiFileAchVecs, SimilarityMat, TrainingSets, TestingSets)

optimal_tau_t = MultiFileAchVecs{1}(~cellfun(@isempty,{MultiFileAchVecs{1}.is_optimal_tau}));
fileInfo = optimal_tau_t(([MultiFileAchVecs{1}.is_optimal_tau] == 1)).dataInfo.FileInfo;

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
            
        ClustersMat = get_clusters_mat(ClusteringDataSets{iCrossValid}(iTau), SimilarityMat(iTau));
        TrainingClusterSets{iCrossValid}(iTau) = prepare_clusters(MultiFileAchVecs, SimilarityMat, ClustersMat, ClusteringDataSets{iCrossValid}, TrainingSets{iCrossValid}, fileInfo, iTau, 1);
        TestingClusterSets{iCrossValid}(iTau) = prepare_clusters(MultiFileAchVecs, SimilarityMat, ClustersMat, ClusteringDataSets{iCrossValid}, TestingSets{iCrossValid}, fileInfo, iTau, 0);       
    end
end

close(f_waitbar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetClusters = prepare_clusters(MultiFileAchVecs, SimilarityMat, ClustersMat, ClusteringData, TestingSet, fileInfo, iTau, isTrain)

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
                    find_testing_avalanche_cluster(ach_vectors_t(iVec), ClusteringData(iTau), SimilarityMat(iTau), ClustersMat);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClustersMat = get_clusters_mat(ClusteringData, SimilarityMat)

ClustersMat = [];

nofLen = length(ClusteringData.Clusters);
for iLen = 1:nofLen
    
    clusters = unique(ClusteringData.Clusters{iLen});
    if isempty(clusters)
        continue;
    end
    
    for iClust = clusters'
        if iLen < nofLen
            ClustersMat{iLen}(iClust).origLen = iLen;
            ClustersMat{iLen}(iClust).origClust = iClust;
            ClusterAchIds = ClusteringData.Id{iLen}(iClust == ClusteringData.Clusters{iLen});
            ClustersMat{iLen}(iClust).idx = contains(SimilarityMat.Id{iLen}, ClusterAchIds);
            clustMat = SimilarityMat.Mat{iLen}(ClustersMat{iLen}(iClust).idx,ClustersMat{iLen}(iClust).idx);
            if numel(clustMat) > 1
                ClustersMat{iLen}(iClust).mean = mean(clustMat(~eye(size(clustMat))));
            else
                ClustersMat{iLen}(iClust).mean = [];
            end
        else %concat
            ClustersMat{iLen}(iClust).origLen = floor(iClust/ClusteringData.concat_cluster_id_prefix_Len);
            ClustersMat{iLen}(iClust).origClust = mod(iClust,ClusteringData.concat_cluster_id_prefix_Len);
            ClustersMat{iLen}(iClust).idx = ClustersMat{ClustersMat{iLen}(iClust).origLen}(ClustersMat{iLen}(iClust).origClust).idx;
            ClustersMat{iLen}(iClust).mean = ClustersMat{ClustersMat{iLen}(iClust).origLen}(ClustersMat{iLen}(iClust).origClust).mean;   
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
function [cluster_num, cluster_sim] = find_testing_avalanche_cluster(avalanche_vec, ClusteringData, SimilarityMat, ClustersMat)

nofLen = length(ClusteringData.Clusters);
cluster_num = zeros(1,nofLen);
cluster_sim = zeros(1,nofLen);

lenInx = [avalanche_vec.length_bins  nofLen-1  nofLen];
for iLen = lenInx
    
    clusters = unique(ClusteringData.Clusters{iLen});
    if isempty(clusters)
        continue;
    end
    
    max_sim = zeros(1,max(clusters));
    mean_sim = zeros(1,max(clusters));
    normalized_mean_sim = zeros(1,max(clusters));
    for iClust = clusters'
        if iLen < nofLen
            ach_idx = contains(SimilarityMat.Id{iLen}, avalanche_vec.id);
            similarities = SimilarityMat.Mat{iLen}(ClustersMat{iLen}(iClust).idx,ach_idx);
        else %concat
            if ClustersMat{iLen}(iClust).origLen ~= avalanche_vec.length_bins
                continue;
            end
            ach_idx = contains(SimilarityMat.Id{avalanche_vec.length_bins}, avalanche_vec.id);
            similarities = SimilarityMat.Mat{avalanche_vec.length_bins}(ClustersMat{iLen}(iClust).idx,ach_idx);
        end
        max_sim(iClust) = max(similarities);
        if max_sim(iClust) == 1 %exact match -> cluster found
            mean_sim(iClust) = 1;
            normalized_mean_sim(iClust) = Inf;
            break;
        else
            mean_sim(iClust) = mean(similarities);
            if ~isempty(ClustersMat{iLen}(iClust).mean) && ClustersMat{iLen}(iClust).mean > 0
                normalized_mean_sim(iClust) = mean_sim(iClust)/ClustersMat{iLen}(iClust).mean;
            else
                normalized_mean_sim(iClust) = mean_sim(iClust)/mean([ClustersMat{iLen}.mean]);
                if isnan(normalized_mean_sim(iClust)) || isinf(normalized_mean_sim(iClust))
                    normalized_mean_sim(iClust) = mean_sim(iClust);
                end
            end
        end
    end
    [cluster_sim(iLen),cluster_num(iLen)] = max(normalized_mean_sim);%max(mean_sim);max(max_sim);
    cluster_sim(iLen) = mean_sim(cluster_num(iLen)); %max_sim(cluster_num(iLen));

end
