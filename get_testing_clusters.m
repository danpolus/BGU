%
%Step 6
% find clusters for testing set vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_testing_clusters()

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';

params_t = global_params();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fn, fp] = uigetfile([fp '*.mat'], 'Select clustering results file');


load([fp fn],'ClusteringData','MultiFileAchVecs','SimilarityMat','TestingSet');

tic

TestingSetClusters = [];
for iTau = 1:length(ClusteringData)
    if isempty(ClusteringData(iTau).tau)
        continue;
    end
    
    TestingSetClusters(iTau).tau = TestingSet(iTau).tau;
    TestingSetClusters(iTau).is_optimal_tau = TestingSet(iTau).is_optimal_tau;
    TestingSetClusters(iTau).CondIds = TestingSet(iTau).CondIds;
    
    for iMode = 1:length(params_t.compare_modes)
        for iCond = 1:length(TestingSet(iTau).CondIds)
            %test epochs vectors
            for iEpoch = 1:length(TestingSet(iTau).EpochIds{iCond})
                %fln = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(4:6));
                %epc = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(10:12));
                ach_vectors_t = MultiFileAchVecs{str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(4:6))}(iTau).epochs_vecs{str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(10:12))};
                
                for iVec = 1:length(ach_vectors_t)
                    avch_length_bins = length(ach_vectors_t(iVec).vec)/MultiFileAchVecs{1}(iTau).nof_channels;
                    [cluster_num, cluster_sim] = find_vector_cluster(ach_vectors_t(iVec).vec, avch_length_bins, MultiFileAchVecs, SimilarityMat, ClusteringData, params_t.similarity_method, params_t.compare_modes{iMode}, iTau, 0);
                    TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_num(iVec) = cluster_num;
                    TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).cluster_sim(iVec) = cluster_sim;
                    TestingSetClusters(iTau).(params_t.compare_modes{iMode})(iCond).EpochClst(iEpoch).avch_length_bins(iVec) = avch_length_bins;
                end
            end
        end
    end
end

toc

save([fp fn(1:end-4) '_testclust.mat'],'TestingSetClusters','ClusteringData');


                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cluster_num, cluster_sim] = find_vector_cluster(v, avch_length_bins, MultiFileAchVecs, SimilarityMat, ClusteringData, similarity_method, compare_mode, iTau, plot_flg)

switch compare_mode
    case 'Full'
        clusters_v = ClusteringData(iTau).ClustersFull;
        Id_list = SimilarityMat(iTau).Id;
    case 'Len'
        clusters_v = ClusteringData(iTau).ClustersLen{avch_length_bins};
        Id_list = SimilarityMat(iTau).IdLen{avch_length_bins};
    case 'ConcatLen'
        clusters_v = ClusteringData(iTau).ClustersConcatLen;
        Id_list = ClusteringData(iTau).IdLConcatLen;
end

clusters = unique(clusters_v);
clusters_similarity = zeros(1,length(clusters));
if plot_flg && ~strcmp(compare_mode, 'Full')
    plot_hist = zeros(length(clusters),length(v));
end
for iClust = 1:length(clusters)
    ClusterAchIds = Id_list(clusters(iClust) == clusters_v);
    if strcmp(compare_mode, 'ConcatLen')
        v_compare = MultiFileAchVecs{str2num(ClusterAchIds{1}(21:23))}(iTau).epochs_vecs{str2num(ClusterAchIds{1}(27:29))}(str2num(ClusterAchIds{1}(33:36))).vec;
        if  avch_length_bins ~= length(v_compare)/MultiFileAchVecs{1}(iTau).nof_channels
            continue;
        end
    end
    for iAch = 1:length(ClusterAchIds)
        %fln = str2num(ClusterAchIds{iAch}(21:23));
        %epc = str2num(ClusterAchIds{iAch}(27:29));
        %ach = str2num(ClusterAchIds{iAch}(33:36));
        v_compare = MultiFileAchVecs{str2num(ClusterAchIds{iAch}(21:23))}(iTau).epochs_vecs{str2num(ClusterAchIds{iAch}(27:29))}(str2num(ClusterAchIds{iAch}(33:36))).vec;
        if strcmp(compare_mode, 'Full')
            %clusters_similarity(iClust) = clusters_similarity(iClust) + calc_sliding_similarity(v, v_compare, MultiFileAchVecs{1}(iTau).nof_channels, similarity_method); %average
            clusters_similarity(iClust) = max(clusters_similarity(iClust), calc_sliding_similarity(v, v_compare, MultiFileAchVecs{1}(iTau).nof_channels, similarity_method)); %nearest neighbour
        else
            %clusters_similarity(iClust) = clusters_similarity(iClust) + calc_vectors_similarity(v, v_compare, similarity_method); %average
            clusters_similarity(iClust) = max(clusters_similarity(iClust), calc_vectors_similarity(v, v_compare, similarity_method)); %nearest neighbour
            if plot_flg
                plot_hist(iClust,:) = plot_hist(iClust,:) + v_compare;
            end
        end
        if clusters_similarity(iClust) == 1 % >0.7 %nearest neighbour
            %figure;stem([v_compare; v]');legend('nearest cluster vector','vector');
            break;
        end
    end
    %clusters_similarity(iClust) = clusters_similarity(iClust)/length(ClusterAchIds); %average 
    if clusters_similarity(iClust) == 1
        break;
    end
end
[cluster_sim,max_inx] = max(clusters_similarity);
cluster_num = clusters(max_inx);

if plot_flg && ~strcmp(compare_mode, 'Full')
    figure;stem([plot_hist(max_inx,:)/max(plot_hist(max_inx,:))/0.5; v]');legend('cluster histogram','vector');title(['cluster: ' num2str(cluster_num) '   similarity: ' num2str(cluster_sim)]);
end
