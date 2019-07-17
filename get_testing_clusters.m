%
%Step 6 : get_testing_clusters - find clusters for testing set vectors
%
%  inputs:
% ClusteringData - clusters and statistics
% MultiFileAchVecs - avalanche vectors from all sets
% TestingSet - avalanches (epochs) for testing
% saveFlg
%
%  outputs:
% TestingSetClusters - testing set clusters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TestingSetClusters = get_testing_clusters(ClusteringData, MultiFileAchVecs, TestingSet, saveFlg)

params_t = global_params();

fileInfo = MultiFileAchVecs{1}(~cellfun(@isempty,{MultiFileAchVecs{1}.is_optimal_tau})).dataInfo.FileInfo;

if saveFlg
    output_fp = [fileInfo.base_fp '5 testing\'];
    mkdir(output_fp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

TestingSetClusters = [];
for iTau = 1:length(ClusteringData)
    
    if isempty(ClusteringData(iTau).tau)
        continue;
    end
    
    TestingSetClusters(iTau).CondClst = [];
    TestingSetClusters(iTau).CondIds = TestingSet(iTau).CondIds;
    TestingSetClusters(iTau).tau = TestingSet(iTau).tau;
    TestingSetClusters(iTau).is_optimal_tau = TestingSet(iTau).is_optimal_tau;
    TestingSetClusters(iTau).fileInfo.base_fp = fileInfo.base_fp;
    TestingSetClusters(iTau).fileInfo.orig_fn = fileInfo.orig_fn;
    
    for iCond = 1:length(TestingSet(iTau).CondIds)
        %test epochs vectors
        for iEpoch = 1:length(TestingSet(iTau).EpochIds{iCond})
            %fln = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(4:6));
            %epc = str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(10:12));
            ach_vectors_t = MultiFileAchVecs{str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(4:6))}(iTau).epochs_vecs{str2num(TestingSet(iTau).EpochIds{iCond}{iEpoch}(10:12))};
            for iVec = 1:length(ach_vectors_t)
                [TestingSetClusters(iTau).CondClst(iCond).EpochClst(iEpoch).cluster_num(iVec,:), TestingSetClusters(iTau).CondClst(iCond).EpochClst(iEpoch).cluster_sim(iVec,:)] = ...
                    find_vector_cluster(ach_vectors_t(iVec).vec, ach_vectors_t(iVec).length_bins, MultiFileAchVecs, ClusteringData, params_t.similarity_method, iTau, 0);
            end
        end
    end
    
end

display(['get testing clusters: ' fileInfo.orig_fn]);
toc

if saveFlg
    save([output_fp fileInfo.orig_fn '_testclust.mat'],'TestingSetClusters','ClusteringData');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cluster_num, cluster_sim] = find_vector_cluster(v, avch_length_bins, MultiFileAchVecs, ClusteringData, similarity_method, iTau, plotFlg)

nofLen = length(ClusteringData(iTau).Clusters);
cluster_num = zeros(1,nofLen);
cluster_sim = zeros(1,nofLen);

lenInx = [avch_length_bins  nofLen-1  nofLen];
for iLen = lenInx
    
    clusters = unique(ClusteringData(iTau).Clusters{iLen});
    if isempty(clusters)
        continue;
    end
    clusters_similarity = zeros(1,length(clusters));
    
    if plotFlg && iLen ~= nofLen-1 %not full matrix
        plot_hist = zeros(length(clusters),length(v));
    end
    
    for iClust = 1:length(clusters)
        ClusterAchIds = ClusteringData(iTau).Id{iLen}(clusters(iClust) == ClusteringData(iTau).Clusters{iLen});
        if iLen == nofLen %concat
            v_compare = MultiFileAchVecs{str2num(ClusterAchIds{1}(21:23))}(iTau).epochs_vecs{str2num(ClusterAchIds{1}(27:29))}(str2num(ClusterAchIds{1}(33:36))).vec;
            if  avch_length_bins ~= length(v_compare)/MultiFileAchVecs{1}(iTau).nof_channels
                continue;
            end
        end
        
        %fln = str2num(ClusterAchIds{iAch}(21:23));
        %epc = str2num(ClusterAchIds{iAch}(27:29));
        %ach = str2num(ClusterAchIds{iAch}(33:36));

        if strcmp(similarity_method, 'jaccard') && iLen ~= nofLen-1 %not full matrix : faster computation
            Vcompare = [];
            for iAch = 1:length(ClusterAchIds)
                Vcompare = [Vcompare;  MultiFileAchVecs{str2num(ClusterAchIds{iAch}(21:23))}(iTau).epochs_vecs{str2num(ClusterAchIds{iAch}(27:29))}(str2num(ClusterAchIds{iAch}(33:36))).vec];
            end            
            sim_M = squareform(1 - pdist([v; Vcompare],'jaccard'), 'tomatrix');
            %clusters_similarity(iClust) = sum(sim_M(1,:)); %average
            clusters_similarity(iClust) = max(sim_M(1,:)); %nearest neighbour
        else
            for iAch = 1:length(ClusterAchIds)
                v_compare = MultiFileAchVecs{str2num(ClusterAchIds{iAch}(21:23))}(iTau).epochs_vecs{str2num(ClusterAchIds{iAch}(27:29))}(str2num(ClusterAchIds{iAch}(33:36))).vec;
                if iLen == nofLen-1 %full
                    %clusters_similarity(iClust) = clusters_similarity(iClust) + calc_sliding_similarity(v, v_compare, MultiFileAchVecs{1}(iTau).nof_channels, similarity_method); %average
                    clusters_similarity(iClust) = max(clusters_similarity(iClust), calc_sliding_similarity(v, v_compare, MultiFileAchVecs{1}(iTau).nof_channels, similarity_method)); %nearest neighbour
                else
                    %clusters_similarity(iClust) = clusters_similarity(iClust) + calc_vectors_similarity(v, v_compare, similarity_method); %average
                    clusters_similarity(iClust) = max(clusters_similarity(iClust), calc_vectors_similarity(v, v_compare, similarity_method)); %nearest neighbour
                    
                    if plotFlg
                        plot_hist(iClust,:) = plot_hist(iClust,:) + v_compare;
                    end
                end
                if clusters_similarity(iClust) == 1 % >0.7 %nearest neighbour
                    %figure;stem([v_compare; v]');legend('nearest cluster vector','vector');
                    break;
                end
            end
        end
        %clusters_similarity(iClust) = clusters_similarity(iClust)/length(ClusterAchIds); %average
        
        if clusters_similarity(iClust) == 1 %for iClust
            break;
        end
    end %for iClust
    
    [cluster_sim(iLen),max_inx] = max(clusters_similarity);
    cluster_num(iLen) = clusters(max_inx);
    
    if plotFlg && iLen ~= nofLen-1 %not full matrix
        figure;stem([plot_hist(max_inx,:)/max(plot_hist(max_inx,:))/0.5; v]');legend('cluster histogram','vector');title(['cluster: ' num2str(cluster_num(iLen)) '   similarity: ' num2str(cluster_sim(iLen))]);
    end
    
end
